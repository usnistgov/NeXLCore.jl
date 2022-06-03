# Define the default algorithms for various pieces of X-ray physics data.
# If you want to change the default data that NeXL uses, make the change in here.

import Unitful: @u_str, ustrip
import PhysicalConstants.CODATA2018: PlanckConstant, SpeedOfLightInVacuum

const plancksConstant = ustrip(PlanckConstant |> u"eV*s")
const hc = ustrip((PlanckConstant * SpeedOfLightInVacuum) |> u"eV*cm") # (plancks⋅speed-of-light)
const speedOfLight = ustrip(SpeedOfLightInVacuum |> u"cm/s")



"""
     energy(ass::AtomicSubShell)
     energy(elm::Element, ss::SubShell)

 The edge energy in eV for the specified AtomicSubShell

Example:

    julia> energy(n"Fe L3")
    708.0999999999999
    julia> energy(n"Fe", n"L3")
    708.0999999999999

    energy(elm::Element, tr::Transition)::Float64
    energy(cxr::CharXRay)

The characteristic X-ray energy for the specified element / transition or characteristic X-ray.
"""
energy(ass::AtomicSubShell)::Float64 = edgeenergy(ass.z, ass.subshell.index)
energy(elm::Element, ss::SubShell)::Float64 = edgeenergy(z(elm), ss.index, alg)
energy(elm::Element, tr::Transition)::Float64 = xrayenergy(z(elm), tr.innershell.index, tr.outershell.index)
energy(cxr::CharXRay) = xrayenergy(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index)

"""
    eachelement()

Return the range of atomic numbers for which there is a complete set of energy, weight, MAC, ... data
"""
eachelement

let allelements = NTuple{length(99), Element}( elements[1:99] )
    global eachelement() = allelements
end

"""
    mac(elm::Element, energy::Float64)::Float64
    mac(elm::Element, cxr::CharXRay)::Float64

The mass absorption coefficient for an X-ray of the specified energy (eV) or
characteristic X-ray line in the specified element.  In cm²/g.
"""
mac(elm::Element, energy::Float64)::Float64 = mac(z(elm), energy)
mac(elm::Element, cxr::CharXRay)::Float64 = mac(z(elm), z(cxr), innerindex(cxr), outerindex(cxr))

"""
    macU(elm::Element, energy::Float64)
    macU(elm::Element, cxr::Float64)::UncertainValue
    macU(elm::Element, cxr::CharXRay)::UncertainValue

The mass absorption coefficient (with uncertainty estimate) for an X-ray of the specified energy (eV) 
or characteristix X-ray line in the specified element.
"""
macU(elm::Element, cxr::CharXRay)::UncertainValue = macU(elm, energy(cxr), alg)

"""
    characteristicXRayAvailable(z::Int, inner::Int, outer::Int)::Float64

Is the weight associated with this transition greater than zero?
"""
charactericXRayAvailable(z::Int, inner::Int, outer::Int)::Bool =
    xrayweight(NormalizeRaw, z, inner, inner, outer) > 0.0

"""
    ionizationcrosssection(ass::AtomicSubShell, energy::AbstractFloat, ty::Type{<:NeXLAlgorithm}=Bote2009)
    ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat)

Computes the absolute ionization crosssection (in cm²) for the specified AtomicSubShell and
electon energy (in eV) using the default algorithm.

Example:

    julia> (/)(map(e->NeXLCore.ionizationcrosssection(n"Fe K",e),[10.0e3,20.0e3])...)
    0.5672910174711278
"""
ionizationcrosssection(
    ass::AtomicSubShell,
    energy::AbstractFloat,
    ty::Type{<:NeXLAlgorithm} = Bote2009,
) = ionizationcrosssection(ass.z, ass.subshell.index, energy, ty)
ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat) =
    ionizationcrosssection(z, ss, energy, Bote2009)
