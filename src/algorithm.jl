# Define the default algorithms for various pieces of X-ray physics data.
# If you want to change the default data that NeXL uses, make the change in here.

import Unitful: @u_str, ustrip
import PhysicalConstants.CODATA2018: PlanckConstant, SpeedOfLightInVacuum

const plancksConstant = ustrip(PlanckConstant |> u"eV*s")
const hc = ustrip((PlanckConstant * SpeedOfLightInVacuum) |> u"eV*cm") # (plancks⋅speed-of-light)
const speedOfLight = ustrip(SpeedOfLightInVacuum |> u"cm/s")


"""
    edgeenergy(z::Int, ss::Int, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64
    edgeenergy(cxr::CharXRay, ::Type{<:NeXLAlgorithm}=FFASTDB)

Return the minimum energy (in eV) necessary to ionize the specified sub-shell in the specified atom
or the ionized shell for the specified characteristic X-ray.
"""
edgeenergy(z::Int, ss::Int)::Float64 = 
    edgeenergy(z, ss, FFASTDB)
edgeenergy(cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB) =
    edgeenergy(cxr.z, cxr.transition.innershell.index, alg)

hasedge(z::Int, ss::Int) = hasedge(z, ss, FFASTDB)

"""
     energy(ass::AtomicSubShell, ty::Type{<:NeXLAlgorithm}=FFASTDB)
     energy(elm::Element, ss::SubShell, ty::Type{<:NeXLAlgorithm}=FFASTDB)

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
energy(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    edgeenergy(ass.z, ass.subshell.index, alg)
energy(elm::Element, ss::SubShell, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    edgeenergy(z(elm), ss.index, alg)
energy(elm::Element, tr::Transition, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64 =
    energy(z(elm), tr.innershell.index, tr.outershell.index, alg)
energy(cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB) = energy(
        cxr.z, 
        cxr.transition.innershell.index,
        cxr.transition.outershell.index,
        alg
    )
    

"""
    eachelement(alg::Type{<:NeXLAlgorithm} = FFASTDB)

Return the range of atomic numbers for which there is a complete set of energy, weight, MAC, ... data
"""
eachelement() = eachelement(FFASTDB)

"""
    strength(cxr::CharXRay)::Float64

The fraction of ionizations of `inner(cxr)` that relax via a characteristic X-ray resulting
from an electronic transition from `outer(cxr)` to `inner(cxr)`.

See also `weight(cxr)`.
"""
strength(cxr::CharXRay)::Float64 = strength(element(cxr), cxr.transition, CullenEADL)


"""
    subshellindices(z::Int, alg::Type{<:NeXLAlgorithm} = FFASTDB)

Return the shells occupied in a neutral, ground state atom of the specified atomic number.
"""
subshellindices(z::Int) = subshellindices(z, FFASTDB)


"""
    jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) =

Compute the jump ratio.
"""
jumpratio(z::Int, ss::Int) = jumpratio(z, ss, FFAST)

"""
    mac(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64
    mac(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::Float64

The mass absorption coefficient for an X-ray of the specified energy (eV) or
characteristic X-ray line in the specified element.  In cm²/g.
"""
mac(elm::Element, energy::Float64)::Float64 =
    mac(elm, energy, FFASTDB)

const userMacs = Dict{Tuple{Element, CharXRay}, Float64}()
mac(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm}=FFASTDB)::Float64 =
    get(userMacs, (elm, cxr), mac(elm, energy(cxr), alg))
set_user_mac!(elm::Element, cxr::CharXRay, mac::Float64) = 
    userMacs[(elm, cxr)]=mac
delete_user_mac!(elm::Element, cxr::CharXRay) = 
    delete!(userMacs, (elm, cxr))
clear_user_macs!() = 
    empty!(userMacs)

"""
    macU(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)
    macU(elm::Element, cxr::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue
    macU(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue

The mass absorption coefficient (with uncertainty estimate) for an X-ray of the specified energy (eV) 
or characteristix X-ray line in the specified element.
"""
macU(elm::Element, energy::Float64)::UncertainValue =
    macU(elm, energy, FFASTDB)
macU(elm::Element, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB)::UncertainValue =
    macU(elm, energy(cxr), alg)

"""
    fluorescenceyield(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

The fraction of `inner` sub-shell ionizations that relax via a characteristic X-ray resulting from an
electronic transition from `outer` to `inner`.
"""
fluorescenceyield(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64 =
    totalWeight(z, inner, inner, outer, alg)

"""
    characteristicyield(z::Int, ionized::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

The fraction of `ionized` sub-shell ionizations that relax via a characteristic X-ray resulting from an
electronic transition from `outer` to `inner`.  This includes both direct transitions (where `outer`==`ionized`)
and cascade (where `outer` != `ionized` due to Coster-Kronig and previous decays.)
"""
characteristicyield(z::Int, ionized::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64 =
    totalWeight(z, ionized, inner, outer, alg)

"""
    characteristicXRayAvailable(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

Is the weight associated with this transition greater than zero?
"""
charactericXRayAvailable(z::Int, inner::Int, outer::Int, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Bool =
    isAvailable(z, inner, outer, alg)

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

    """
    strength(elm::Element, tr::Transition)::Float64

Return the nominal line strenth for the specified transition in the specified element.
The strength differs from the weight by the fluorescence yield.  Assumes an overvoltage of 4.0
"""
strength(elm::Element, tr::Transition, ty::Type{<:NeXLAlgorithm} = CullenEADL)::Float64 =
    ionizationfraction(z(elm), tr.innershell.index, 4.0) *
    fluorescenceyield(z(elm), tr.innershell.index, tr.outershell.index, ty)


"""
    fluorescenceyield(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.  Does not include Coster-Kronig
"""
fluorescenceyield(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64 = sum(map(
    s -> fluorescenceyield(ass.z, ass.subshell.index, s, alg),
    ass.subshell.index+1:length(allsubshells),
))

"""
    fluorescenceyieldcc(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64

The fraction of relaxations from the specified shell that decay via radiative transition
rather than electronic (Auger) transition.  Includes Coster-Kronig
"""
function fluorescenceyieldcc(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64
    f(ss) = sum(map(
        s -> fluorescenceyield(ass.z, ass.subshell.index, s, alg),
        ss.index+1:length(allsubshells),
    ))
    return sum(map(ss -> f(ss), ass.subshell.index+1:lastsubshell(shell(ass)).index))
end
