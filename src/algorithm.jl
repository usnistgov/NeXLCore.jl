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
    

Base.:(==)(elm1::Element, elm2::Element) = z(elm1) == z(elm2)

"""
    eachelement()

Return the range of atomic numbers for which there is a complete set of energy, weight, MAC, ... data
"""
eachelement

let allelements = NTuple{length(FFAST.eachelement()), Element}( elements[FFAST.eachelement()] )
    global eachelement() = allelements
end

"""
    subshellindices(z::Int, alg::Type{<:NeXLAlgorithm} = FFASTDB)

Return the shells occupied in a neutral, ground state atom of the specified atomic number.
"""
subshellindices(z::Int) = subshellindices(z, FFASTDB)


"""
    jumpratio(z::Int, ss::Int) =

Compute the jump ratio.
"""
jumpratio(z::Int, ss::Int) = jumpratio(z::Int, ss::Int, FFASTDB)

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
Represents the fractional number of X-rays emitted following the ionization of the sub-shell `ionized` via
the characteristic X-ray `z inner-outer`.  Due to cascades, `inner` does not necessarily equal `ionized`.
The `ionized` subshell may transition to a valency in `inner` via a combination of Auger, fluorescence or
Koster-Kronig transitions.  The various different forms make assumptions about the relationship between
`ionized` and `inner`, and about `outer`.

    fluorescenceyield(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL)::Float64

The fraction of relaxations from the specified shell that relax via any radiative transition. (`inner`==`ionized`)

    fluorescenceyield(cxr::CharXRay, alg::Type{<:NeXLAlgorithm}=CullenEADL)

The fraction of ionizations of `inner(cxr)` that relax via the one path `cxr`. `ionized==inner` && outer(cxr)

    fluorescenceyield(ash::AtomicSubShell, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = CullenEADL)::Float64

The fractional number of `cxr` X-rays emitted (on average) for each ionization of `ash`.  This makes no 
assumptions about `inner`, `outer` and `ionized`
"""
fluorescenceyield(ass::AtomicSubShell, alg::Type{<:NeXLAlgorithm}=CullenEADL) =
    fluorescenceyield(ass.z, ass.subshell.index, alg)
fluorescenceyield(cxr::CharXRay, alg::Type{<:NeXLAlgorithm}=CullenEADL) =
    fluorescenceyield(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index, alg)
fluorescenceyield(ash::AtomicSubShell, cxr::CharXRay, alg::Type{<:NeXLAlgorithm} = CullenEADL) =
    ash.z == cxr.z ? fluorescenceyield(ash.z, ash.subshell.index, cxr.inner.index, cxr.outer.index, alg) : 0.0
