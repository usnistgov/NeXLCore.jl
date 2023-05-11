import Unitful: @u_str, ustrip
import PhysicalConstants.CODATA2018: PlanckConstant, SpeedOfLightInVacuum
import BoteSalvatICX # For ionization crosssections

# These evaluate at compile time...
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
energy(ass::AtomicSubShell, ::Type{DefaultAlgorithm})::Float64 = edgeenergy(ass.z, ass.subshell.index)
energy(elm::Element, ss::SubShell, ::Type{DefaultAlgorithm})::Float64 = edgeenergy(z(elm), ss.index)
energy(elm::Element, tr::Transition, ::Type{DefaultAlgorithm})::Float64 = xrayenergy(z(elm), tr.innershell.index, tr.outershell.index)
energy(cxr::CharXRay, ::Type{DefaultAlgorithm}) = xrayenergy(cxr.z, cxr.transition.innershell.index, cxr.transition.outershell.index)

energy(ass::AtomicSubShell)::Float64 = energy(ass, DefaultAlgorithm)
energy(elm::Element, ss::SubShell)::Float64 = energy(elm, ss, DefaultAlgorithm)
energy(elm::Element, tr::Transition)::Float64 = energy(elm, tr, DefaultAlgorithm)
energy(cxr::CharXRay) = energy(cxr::CharXRay, DefaultAlgorithm)

"""
    edgeenergy(cxr::CharXRay)

Returns the energy associated with the inner shell of this characteristic X-ray (eV).
"""
edgeenergy(cxr::CharXRay) = edgeenergy(z(cxr), innerindex(cxr))

"""
    eachelement()

Return the range of atomic numbers for which there is a complete set of energy, weight, MAC, ... data
"""
global eachelement() = elements[allz]

"""
    mac(elm::Element, energy::Float64)::Float64
    mac(elm::Element, cxr::CharXRay)::Float64

The mass absorption coefficient for an X-ray of the specified energy (eV) or
characteristic X-ray line in the specified element.  In cm²/g.
"""
mac(elm::Element, energy::Float64, ::Type{DefaultAlgorithm})::Float64 = mac(z(elm), energy)
mac(elm::Element, cxr::CharXRay, ::Type{DefaultAlgorithm})::Float64 = mac(z(elm), z(cxr), innerindex(cxr), outerindex(cxr))
mac(elm::Element, cxr::CharXRay, ::Type{DTSA}) = mac(elm, energy(cxr), DTSA)

mac(elm::Element, energy::Float64) = mac(elm, energy, DefaultAlgorithm)
mac(elm::Element, cxr::CharXRay) = mac(elm, cxr, DefaultAlgorithm)

"""
    setmac!(elm::Element, cxr::CharXRay, mac::Float64)

Specify a custom mass absorption coefficient (MAC) for the specified X-ray in the specified element.
"""
setmac!(elm::Element, cxr::CharXRay, mac::Float64) = setmac!(z(elm), z(cxr), innerindex(cxr), outerindex(cxr), mac)

"""
    ressetmac!(elm::Element, cxr::CharXRay)

Restore the default mass absorption coefficient for the specified element and characteristic X-ray.
"""
resetmac!(elm::Element, cxr::CharXRay) = resetmac!(z(elm), z(cxr), innerindex(cxr), outerindex(cxr))

"""
    loadcustommacs!(source::AbstractString, elms::AbstractArray{Element})
    loadcustommacs!(source::AbstractString, elm::Element)

Load custom macs from the database associated with the specified source.

Sources include "Henke1974", "Henke1982", "Bastin1989", "Henke1993", "Bastin1997", "Ruste1979", "Kohlhaas1970", "Weisweiler1975", 
"Bastin1990", "Bastin1988", "Poml2020", "Ruste1975", "Henke1982", "Farthing1990", "Sabbatucci2016"

Use `listcustommacs(cxr)` to explore available MACs.
"""
loadcustommacs!(source::AbstractString, elms::AbstractArray{Element}) = loadcustommacs!(source, z.(elms))
loadcustommacs!(source::AbstractString, elm::Element) = loadcustommacs!(source, [ z(elm) ])

"""
    loadcustommac!(elm::Element, cxr::CharXRay, source::AbstractString) 

Load the custom mass absorption coefficient associated with the specified `CharXRay` in the `Element` from the source.
Use `listcustommacs(cxr)` to explore available MACs.
"""
function loadcustommac!(elm::Element, cxr::CharXRay, source::AbstractString) 
    try
        loadcustommac!(z(elm), z(cxr), innerindex(cxr), outerindex(cxr), source)
    catch
        error("$source does not provide a custom MAC for $cxr in $(symbol(elm)).")
    end
end


"""
    listcustommacs(cxr::CharXRay)
Generated a list of available custom MACs for `cxr` in various elements.

    listcustommacs(elms::Set{Element}|AbstractVector{Element}|Element...|Material)

Generate a list of available custom MACS for elements and X-rays both produced and absorbed by these elements.
"""
function listcustommacs(cxr::CharXRay)
    res = listcustommacs(z(cxr), innerindex(cxr), outerindex(cxr))
    sort!(res, lt=(a,b)->a[2]<b[2])
    map(res) do (ref, z1, z2, inner, outer, MACpi) 
        cxr2 = CharXRay(z2, Transition(SubShell(inner), SubShell(outer)))
        "MAC[$cxr2 in $(symbol(elements[z1])), $ref] = $MACpi"
    end
end

function listcustommacs(elms::AbstractSet{Element})
    inmat(r) = elements[r[3]] in elms
    mapreduce(append!, elms) do elm
        map(filter(inmat, listcustommacs(z(elm)))) do (ref, z1, z2, inner, outer, MACpi)
            cxr2 = CharXRay(z2, Transition(SubShell(inner), SubShell(outer)))
            "MAC[$cxr2 in $(symbol(elements[z1])), $ref] = $MACpi"
        end
    end
end
listcustommacs(elms::AbstractVector{Element}) = listcustommacs(Set(elms))
listcustommacs(elms::Element...) = listcustommacs(Set(elms))

abstract type MACUncertainty end
struct MonatomicGas <: MACUncertainty end
struct SolidLiquid <: MACUncertainty end

"""
    fractionaluncertainty(::Type{MonatomicGas}, z::Integer, energy)
Determines from the element and energy, the approximate range of fractional uncertainties to
associate with the total and photoelectric components of the mass attenuation coefficients
for monatomic gas samples.
Based on [this table](https://physics.nist.gov/PhysRefData/FFast/Text2000/sec06.html#tab2).
"""
function fractionaluncertainty(::Type{MonatomicGas}, z::Integer, energy)
    low, high = 0.01, 0.01
    if energy < 200.0
        low, high = 0.5, 1.0
    elseif energy < 500.0
        low, high = 0.20, 0.30
    elseif energy < 1.0
        low, high = 0.03, 0.10
    end
    distance = ( (energy - edgeenergy(z, sh)) / energy for sh in eachedge(z) )
    if minimum(abs.(distance)) < 0.001
        low, high = max(low, 0.2), max(high, 0.3)
    end
    u = [ energy / edgeenergy(z, sh) for sh in eachedge(z) ]
    if (u[1] > 1.0) && (u[1] < 1.1)
        low, high = max(low, 0.1), max(high, 0.1)
    elseif (u[1] >= 1.1) && (u[1] < 1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    # L1, M1, M2, M3
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 2, 5, 6, 7 ])
        if u[sh] < 1.15
            low, high = max(low, 0.15), max(high, 0.15)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    # L2, L3, M4, M5
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 3, 4, 8, 9 ])
        if u[sh] < 1.15
            low, high = max(low, 0.20), max(high, 0.20)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy > 200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end

"""
    fractionaluncertainty(::Type{SolidLiquid}, z::Integer, energy)

Determines from the element and energy, the approximate range of fractional uncertainties to
associate with the total and photoelectric components of the mass attenuation coefficients
for solids and liquids.
Based on [this table](https://physics.nist.gov/PhysRefData/FFast/Text2000/sec06.html#tab2).
"""
function fractionaluncertainty(::Type{SolidLiquid}, z::Integer, energy)
    low, high = 0.01, 0.01
    if energy < 200.0
        low, high = 1.0, 2.0
    elseif energy < 500.0
        low, high = 0.50, 1.0
    elseif energy < 1.0
        low, high = 0.05, 0.20
    end
    distance = ( (energy - edgeenergy(z, sh)) / energy for sh in eachedge(z) )
    if minimum(abs.(distance)) < 0.001
        low, high = max(low, 0.5), max(high, 0.5)
    end
    u = [ energy / edgeenergy(z, sh) for sh in eachedge(z) ]
    if (u[1] > 1.0) && (u[1] < 1.1)
        low, high = max(low, 0.1), max(high, 0.2)
    elseif (u[1] >= 1.1) && (u[1] < 1.2)
        low, high = max(low, 0.03), max(high, 0.03)
    end
    # L1, M1, M2, M3
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 2, 5, 6, 7 ])
        if u[sh] < 1.15
            low, high = max(low, 0.15), max(high, 0.30)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    # L2, L3, M4, M5
    for sh in filter(sh->get(u, sh, 0.0) >= 1.0, [ 3, 4, 8, 9 ] )
        if u[sh] < 1.15
            low, high = max(low, 0.20), max(high, 0.40)
        elseif u[sh] < 1.4
            low, high = max(low, 0.04), max(high, 0.04)
        end
    end
    if energy > 200.0e3
        low, high = max(low, 0.02), max(high, 0.03)
    end
    return (low, high)
end

"""
    macU(elm::Element, energy::Float64)
    macU(elm::Element, cxr::CharXRay)::UncertainValue

The mass absorption coefficient (with uncertainty estimate) for an X-ray of the specified energy (eV) 
or characteristix X-ray line in the specified element.
"""
function macU(elm::Element, energy::Float64, alg::Type{<:NeXLAlgorithm}=DefaultAlgorithm, unc::Type{<:MACUncertainty}=SolidLiquid)::UncertainValue
    macv = mac(elm, energy, alg)
    return uv(
        macv,
        min(fractionaluncertainty(unc, z(elm), energy)[1], 0.9) * macv
    )
end
function macU(elm::Element, cxr::CharXRay, unc::Type{<:MACUncertainty}=SolidLiquid)::UncertainValue
    # Force lookup for custom MACs
    macv = mac(z(elm), z(cxr), innerindex(cxr), outerindex(cxr))
    return uv(
        macv,
        min(fractionaluncertainty(unc, z(elm), energy(cxr))[1], 0.9) * macv
    )
end


"""
    characteristicXRayAvailable(z::Int, inner::Int, outer::Int)::Float64

Is the weight associated with this transition greater than zero?
"""
charactericXRayAvailable(z::Int, inner::Int, outer::Int)::Bool = hasxray(z, inner, outer)    

struct Bote2009 <: NeXLAlgorithm end

"""
    ionizationcrosssection(z::Int, shell::Int, energy::AbstractFloat, ::Type{Bote2009})
    ionizationcrosssection(ass::AtomicSubShell, energy::AbstractFloat, ty::Type{<:NeXLAlgorithm}=Bote2009)

Computes the absolute ionization crosssection (in cm²/e⁻) for the specified AtomicSubShell and
electon energy (in eV).

Example:

    julia> (/)(map(e->NeXLCore.ionizationcrosssection(n"Fe K",e),[10.0e3,20.0e3])...)
    0.5672910174711278
"""
function ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat, ::Type{Bote2009})::Float64
    return BoteSalvatICX.hasedge(z, ss) ?
        BoteSalvatICX.ionizationcrosssection(z, ss, energy, BoteSalvatICX.edgeenergy(z, ss)) : #
        0.0
end
ionizationcrosssection(
    ass::AtomicSubShell,
    energy::AbstractFloat,
    ty::Type{<:NeXLAlgorithm} = DefaultAlgorithm
) = ionizationcrosssection(ass.z, ass.subshell.index, energy, ty)
ionizationcrosssection(
    z::Int, #
    ss::Int, #
    energy::AbstractFloat, #
    ty::Type{<:NeXLAlgorithm} = DefaultAlgorithm
) = ionizationcrosssection(z, ss, energy, ty==DefaultAlgorithm ? Bote2009 : ty)


