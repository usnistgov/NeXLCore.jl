"""
η(::Type{<:BackscatterCoefficient}, elm::Element, e0::Real)
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

η(::Type{Tomlin1963}, elm::Element, e0::Real) =
 	log(max(5,z(elm)))/6.0 - 0.25

"""@article{Love_1978,
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

function η(::Type{LoveScott1978η}, elm::Element, e0::Real)
	η20(z) = (-52.3791 + z*(150.48371 + z*(-1.67373 + z*0.00716)))*1.0e-4
	Goη20(z) = (-1112.8 + z*(30.289 + z* -0.15498))*1.0e-4
	zz=convert(Float64, z(elm))
	return η20(zz)*(1.0+Goη20(zz)*log(e0/20.0e3))
end

η(elm::Element, e0::Real) = η(LoveScott1978η, elm, e0)
