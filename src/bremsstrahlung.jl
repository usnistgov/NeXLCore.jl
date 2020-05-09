
using NeXLCore
"""
An abstract structure that implements

    NeXLCore.bremsstrahlung(::Type{<:NeXLBremsstrahlung}, e::AbstractFloat, e0::AbstractFloat, elm::Element; kwargs...)
"""
abstract type NeXLBremsstrahlung <: NeXLAlgorithm end

"""
@article{kramers1923,
    author = { Kramers, H. A. },
    title = {XCIII. On the theory of X-ray absorption and of the continuous X-ray spectrum},
    journal = {The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science},
    volume = {46},
    number = {275},
    pages = {836-871},
    year  = {1923},
    publisher = {Taylor & Francis},
    doi = {10.1080/14786442308565244},
    URL = { https://doi.org/10.1080/14786442308565244 },
    eprint = { https://doi.org/10.1080/14786442308565244}
}
"""
struct Kramers1923 <: NeXLBremsstrahlung end

"""
@InProceedings{ lifshin1974,
    author={Lifshin, E.},
    journal="Proceedings Ninth National Conference on Electron Probe Analysis, Ottowa, CA",
    year={1974},
    volume={174},
    number={3},
}
"""
struct Lifshin1974 <: NeXLBremsstrahlung end

"""
@article{reed1975shape,
  title={The shape of the continuous X-ray spectrum and background corrections for energy-dispersive electron microprobe analysis},
  author={Reed, SJB},
  journal={X-Ray Spectrometry},
  volume={4},
  number={1},
  pages={14--17},
  year={1975},
  publisher={Wiley Online Library}
}
"""
struct Reed1975 <: NeXLBremsstrahlung end

"""
@article{ smithgoldtomlinson1975
    title="The atomic number dependence of the X‐ray continuum intensity and the practical calculation of background in energy dispersive electron microprobe analysis",
    author = {Smith, D. G. W.  and Gold, C. M. and Tomlinson, D. A.},
    journal={X-ray Spectrometry},
    volume={4},
    pages={149-156},
    year={1975},
    doi={10.1002/xrs.1300040311},
}
"""
struct Smith1975 <: NeXLBremsstrahlung end

"""
@article{small1987modeling,
  title={Modeling of the bremsstrahlung radiation produced in pure-element targets by 10--40 keV electrons},
  author={Small, John A and Leigh, Stefan D and Newbury, Dale E and Myklebust, Robert L},
  journal={Journal of applied physics},
  volume={61},
  number={2},
  pages={459--469},
  year={1987},
  publisher={American Institute of Physics}
}
"""
struct Small1987 <: NeXLBremsstrahlung end

"""
@article{trincavelli1998model,
  title={Model for the bremsstrahlung spectrum in EPMA. Application to standardless quantification},
  author={Trincavelli, Jorge and Castellano, Gustavo and Riveros, J Alberto},
  journal={X-Ray Spectrometry: An International Journal},
  volume={27},
  number={2},
  pages={81--86},
  year={1998},
  publisher={Wiley Online Library}
}
"""
struct Trincavelli1997 <: NeXLBremsstrahlung end

"""
@article{castellano2004analytical,
  title={Analytical model for the bremsstrahlung spectrum in the 0.25--20 keV photon energy range},
  author={Castellano, Gustavo and Osan, Janos and Trincavelli, Jorge},
  journal={Spectrochimica Acta Part B: Atomic Spectroscopy},
  volume={59},
  number={3},
  pages={313--319},
  year={2004},
  publisher={Elsevier}
}
"""
struct Castellano2004a <: NeXLBremsstrahlung end
struct Castellano2004b <: NeXLBremsstrahlung end

bremsstrahlung(::Type{Kramers1923}, e::AbstractFloat, e0::AbstractFloat, elm::Element) =
    z(elm)*(e0-e)/e

bremsstrahlung(::Type{Lifshin1974}, e::AbstractFloat, e0::AbstractFloat, elm::Element; a=0.0001) =
    bremsstrahlung(Kramers1923,e,e0,elm)*(1.0-1.0e-3*a*(e0-e))

bremsstrahlung(::Type{Reed1975}, e::AbstractFloat, e0::AbstractFloat, elm::Element; b=0.0001) =
    f*z(elm)*(1.0e-3*(e0-e))/((1.0e-3*e)^(1.0+b))

function bremsstrahlung(::Type{Smith1975}, e::AbstractFloat, e0::AbstractFloat, elm::Element)
    n(e,z) = 1.0e-3*e*(0.0739 - 0.0051*log(z)) + 1.6561 - 0.115*log(z)
    x(e,z) = 1.76 - 0.00145*z/(1.0e-3*e)
    return ((e0-e)/e)^x(e,z(elm))*z^n(e,z(elm))
end

function bremsstrahlung(::Type{Small1987}, e::AbstractFloat, e0::AbstractFloat, elm::Element)
    M(e0) = 0.00599e-3*e0 + 1.05
    B(e0) = -0.0322e-3*e0
    return (z(elm)*((e0-e)/e))^M(e0)*(1.0e-3*e)^B(e0)
end

bremsstrahlung(::Type{Trincavelli1997}, e::AbstractFloat, e0::AbstractFloat, elm::Element) =
    (sqrt(z(elm))*(e0-e)/e)*(-54.86 - 1.072e-3*e + 0.2835e-3*e0 +30.4*log(z(elm)) + 875.0/(z(elm)^2*(1.0e-3*e0)^0.08))

bremsstrahlung(::Type{Castellano2004a}, e::AbstractFloat, e0::AbstractFloat, elm::Element) =
    (sqrt(z(elm))*(e0-e)/e)*(-73.90 - 1.2446e-3*e + 36.502*log(z(elm)) + 148.5e-3*e0^0.1239)* #
            (1.0 + (-0.006624 + 0.002906e-3*e0)*z(elm)/(1.0e-3*e))

function bremsstrahlung(::Type{Castellano2004b}, e::AbstractFloat, e0::AbstractFloat, elm::Element)
    zz, E0, E = z(elm), e0/1000.0, e/1000.0
    return ((-E + E0)*sqrt(zz)*(1.0 + ((-0.006626 + 0.0002906*E0)*zz)/E)*
            (-77.28317356370013 + (148.5*E0^0.1293)/zz + 36.502*log(zz)))/E
end

bremsstrahlung(e::AbstractFloat, e0::AbstractFloat, elm::Element) =
    bremsstrahlung(Castellano2004b, e, e0, elm)

function bremsstrahlung(ty::Type{<:NeXLBremsstrahlung}, e::AbstractFloat, e0::AbstractFloat, mat::Material; kwargs...)
    af=atomicfraction(mat) # According to Trincavelli1997
    return mapreduce(elm->af[elm]*bremmstrahlung(ty,e,e0,elm; kwargs...),+,keys(af))
end
