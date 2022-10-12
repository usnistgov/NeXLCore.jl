"""
    η(::Type{<:BackscatterCoefficient}, elm::Element, e0::Real)

Abstract type to represent backscatter coefficient algorithms.
"""
abstract type BackscatterCoefficient <: NeXLAlgorithm end

"""
@article{Tomlin_1963,
	doi = {10.1088/0370-1328/82/3/118},
	url = {https://doi.org/10.1088%2F0370-1328%2F82%2F3%2F118},
	year = 1963,
	month = {sep},
	publisher = {{IOP} Publishing},
	volume = {82},
	number = {3},
	pages = {465--466},
	author = {S G Tomlin},
	title = {The Back-scattering of Electrons from Solids},
	journal = {Proceedings of the Physical Society},
	abstract = {Archard's diffusion model of electron back-scattering is discussed on the basis of a
	result obtained by Tomlin and Metchnik in 1963 in their treatment of x-ray emission intensities.
	The resulting simple formula for the back-scattering coefficient is in reasonably good agreement
	with measured values.}
}

The implementation adapted to not return negative numbers for z<5.
"""
struct Tomlin1963 <: BackscatterCoefficient end

η(::Type{Tomlin1963}, elm::Element, e0::Float64) = log(max(5, z(elm))) / 6.0 - 0.25

"""
@article{Love_1978,
	doi = {10.1088/0022-3727/11/10/002},
	url = {https://doi.org/10.1088%2F0022-3727%2F11%2F10%2F002},
	year = 1978,
	month = {jul},
	publisher = {{IOP} Publishing},
	volume = {11},
	number = {10},
	pages = {1369--1376},
	author = {G Love and V D Scott},
	title = {Evaluation of a new correction procedure for quantitative electron probe microanalysis},
	journal = {Journal of Physics D: Applied Physics},
	abstract = {A new correction procedure for converting electron-probe microanalysis measurements into
	true weight concentration is proposed. It incorporates a new atomic number correction and an absorption
	correction based upon Bishop's model (1974). Unlike earlier treatments the model does not have to
	rely upon any empirical optimisation of input parameters. The correction procedure has been tested by
	applying it to a wide range of microanalysis data including light-element results, and it is shown to
	give greater accuracy than the established methods.}
}
"""
struct LoveScott1978η <: BackscatterCoefficient end

function η(::Type{LoveScott1978η}, elm::Element, e0::Float64)
    η20(z) = evalpoly(z, (-52.3791e-4, 150.48371e-4, -1.67373e-4, 0.00716e-4))
    Goη20(z) = evalpoly(z, (-1112.8e-4, 30.289e-4,-0.15498e-4))
    zz = Float64(z(elm))
    return η20(zz) * (1.0 + Goη20(zz) * log(e0 / 20.0e3))
end

"""
The model for Pouchou's 1991 model ("Green Book") of the BackscatterCoefficient.
"""
struct Pouchou1991η <: BackscatterCoefficient end

η(::Type{Pouchou1991η}, elm::Element, e0::Float64) =
    0.00175 * z(elm) + 0.37 * (1.0 - exp(-0.015 * z(elm)^1.3))


struct August1989η <: BackscatterCoefficient end

function η(::Type{August1989η}, elm::Element, e0::Float64)
    zz, e0k = Float64(z(elm)), 0.001 * e0
    return evalpoly(log(zz), (0.1904, -0.2236, 0.1292, -0.01491)) *
           (0.9987 + 2.167e-4 * zz) * (e0k^(0.1382 - 0.9211 / sqrt(zz)))
end

"""
The model for Reimer's model of the BackscatterCoefficient.
"""
struct Reimer1998 <: BackscatterCoefficient end

η(::Type{Reimer1998}, elm::Element, e0::Float64) = #
    evalpoly(z(elm), (-0.0254, 0.016, -1.86e-4, 8.3e-7))

η(::Type{Reimer1998}, mat::Material, e0::Float64) = #
    evalpoly(sum(elm->value(mat[elm])*z(elm),keys(mat)), (-0.0254, 0.016, -1.86e-4, 8.3e-7))

"""
    η(::Type{<:BackscatterCoefficient}, mat::Material, e0::Float64) = #
    η(elm::Element, e0::Real)

Models are: Tomlin1963, LoveScott1978η, Pouchou1991η, August1989η, Reimer1998

The default backscatter coefficient algorith is August1989η.
"""
η(elm::Element, e0::Real) = η(August1989η, elm, e0)

"""
    η([ty::Type{<:BackscatterCoefficient},] mat::Material, e0::Float64)::Float64 =

Calculate the backscatter coefficient for a material using Armstrong's 1991 algorithm for materials.

@incollection{armstrong1991quantitative,
  title={Quantitative elemental analysis of individual microparticles with electron beam instruments},
  author={Armstrong, John T},
  booktitle={Electron probe quantitation},
  pages={261--315},
  year={1991},
  publisher={Springer}
}
"""
η(ty::Type{<:BackscatterCoefficient}, mat::Material, e0::Float64)::Float64 =
    mapreduce(el -> η(ty, el, e0) * elasticfraction(el, mat, e0), +, keys(mat))
η(mat::Material, e0::Float64)::Float64 =
    mapreduce(el -> η(el, e0) * elasticfraction(el, mat, e0), +, keys(mat))



"""
    elasticfraction(elm::Element, mat::Material, e0::Float64)::Float64

Computes the fraction of the total scattering cross-section associated with `elm` in `mat` at beam energy `e0`.

@incollection{armstrong1991quantitative,
  title={Quantitative elemental analysis of individual microparticles with electron beam instruments},
  author={Armstrong, John T},
  booktitle={Electron probe quantitation},
  pages={261--315},
  year={1991},
  publisher={Springer}
}
"""
function elasticfraction(elm::Element, mat::Material, e0::Float64)::Float64
    function elasticcrosssection(zz, ek)
        αα, mc2 = 3.4e-3 * zz^0.67 / ek, 511.0
        return 5.21e-21 *
               ((zz / ek)^2) *
               (4π / (αα * (1.0 + αα))) *
               (((ek + mc2) / (ek + 2mc2))^2)
        # return ((511. + ek)*zz^1.33)/((1022. + ek)*(ek + 0.0034*zz^0.67))
        # Which is very close to zz^1.33 for except for large z at low ek
    end
    aw, ek = atomicfraction(mat), 0.001 * e0
    den = mapreduce(el -> value(aw[el]) * elasticcrosssection(z(el), ek), +, keys(mat))
    return value(get(aw, elm, zero(valtype(aw)))) * elasticcrosssection(z(elm), ek) / den
end

"""
	electronfraction(elm::Element, mat::Material)::Float64

The electron fraction as defined in:

@article{donovan2003compositional,
  title={Compositional averaging of backscatter intensities in compounds},
  author={Donovan, John J and Pingitore, Nicholas E and Westphal, Andrew},
  journal={Microscopy and Microanalysis},
  volume={9},
  number={3},
  pages={202--215},
  year={2003},
  publisher={Cambridge University Press}
}
"""
function electronfraction(elm::Element, mat::Material)::Float64
    aw = atomicfraction(mat)
    return value(get(aw, elm, zero(valtype(aw)))) * z(elm) / mapreduce(el -> value(aw[el]) * z(el), +, keys(mat))
end

"""
	zbar(mat::Material)::Float64

The mean atomic number as calculated using

@article{donovan2003compositional,
  title={Compositional averaging of backscatter intensities in compounds},
  author={Donovan, John J and Pingitore, Nicholas E and Westphal, Andrew},
  journal={Microscopy and Microanalysis},
  volume={9},
  number={3},
  pages={202--215},
  year={2003},
  publisher={Cambridge University Press}
}

or, equivalently,

@article{saldick1954backscattering,
  title={Backscattering from Targets of Low Atomic Number Bombarded with 1—2 Mev Electrons},
  author={Saldick, Jerome and Allen, Augustine O},
  journal={The Journal of Chemical Physics},
  volume={22},
  number={10},
  pages={1777--1777},
  year={1954},
  publisher={American Institute of Physics}
}

"""
zbar(mat::Material)::Float64 =
    mapreduce(el -> electronfraction(el, mat) * z(el), +, keys(mat))
