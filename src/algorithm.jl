# Defines how various different algorithms are implemented
import FFAST # for mass absorption coefficienct
import BoteSalvatICX # For ionization crosssections

const plancksConstant = 4.135667696e-15 # eV s
const hc = 1.23984193e-4 # eV cm (plancks⋅speed-of-light)
const speedOfLight = 2.99792458e10 # cm/s

const subshellnames = (
    "K",
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
    "M4",
    "M5",
    "N1",
    "N2",
    "N3",
    "N4",
    "N5",
    "N6",
    "N7",
    "O1",
    "O2",
    "O3",
    "O4",
    "O5",
    "O6",
    "O7",
    "O8",
    "O9",
    "P1",
    "P2",
    "P3",
    "P4",
    "P5",
    "P6",
    "P7",
    "P8",
    "P9",
    "P10",
    "P11",
    "Q1",
    "Q2",
    "Q3",
)
"""
    subshellindex(ss::AbstractString)

Map sub-shell names ("K","L1", ...,"P11") to an integer index ("K"==1,"L1"==2, etc).
"""
subshellindex(ssname::AbstractString) = findfirst(name -> ssname == name, subshellnames)

struct FFASTDB <: NeXLAlgorithm end

"""
    mac(elm::Element, energy::Float64, ::Type{FFASTDB})::Float64

Compute the mass absorption coefficient for the specified energy X-ray (in eV) in the specified element (z=> atomic number).
"""
mac(elm::Element, energy::Float64, ::Type{FFASTDB})::Float64 =
    FFAST.mac(FFAST.PhotoElectricMAC, z(elm), energy)


"""
    macU(elm::Element, energy::Float64, ::Type{FFASTDB})::UncertainValue

Compute the mass absorption coefficient (with estimate of uncertainty) for the specified energy X-ray (in eV) in the specified
element (z=> atomic number).
"""
function macU(elm::Element, energy::Float64, ::Type{FFASTDB})::UncertainValue
    macv = FFAST.mac(FFAST.PhotoElectricMAC, z(elm), energy)
    return uv(
        macv,
        min(FFAST.fractionaluncertainty(FFAST.SolidLiquid, z(elm), energy)[1], 0.9) * macv,
    )
end

"""
    edgeenergy(z::Int, ss::Int, ::Type{FFASTDB})::Float64

Return the minimum energy (in eV) necessary to ionize the specified sub-shell in the specified atom.
"""
edgeenergy(z::Int, ss::Int, ::Type{FFASTDB})::Float64 = FFAST.edgeenergy(z, ss)
edgeenergy(z::Int, ss::Int) = edgeenergy(z, ss, FFASTDB)

"""
    eachelement() = 1:92

Return the range of atomic numbers for which there is a complete set of energy, weight, MAC, ... data
"""
eachelement(::Type{FFASTDB}) = FFAST.eachelement()
eachelement() = eachelement(FFASTDB)
"""
    subshellsindexes(z::Int)

Return the shells occupied in a neutral, ground state atom of the specified atomic number.
"""
subshellsindexes(z::Int, ::Type{FFASTDB}) = FFAST.eachedge(z)
subshellsindexes(z::Int) = subshellsindexes(z::Int, FFASTDB)


"""
    characteristicXRayEnergy(z::Int, inner::Int, outer::Int)::Float64

Return energy (in eV) of the transition by specified inner and outer sub-shell index.
"""
characteristicXRayEnergy(z::Int, inner::Int, outer::Int, ::Type{FFASTDB})::Float64 =
    FFAST.edgeenergy(z, inner) - FFAST.edgeenergy(z, outer)

struct Bote2009 <: NeXLAlgorithm end

"""
    ionizationcrosssection(z::Int, shell::Int, energy::AbstractFloat, ::Type{Bote2009})

Compute the absolute ionization crosssection (in cm²) for the specified element, shell and
electon energy (in eV).
"""
ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat, ::Type{Bote2009}) =
    BoteSalvatICX.hasedge(z, ss) ?
    BoteSalvatICX.ionizationcrosssection(
        z,
        ss,
        energy,
        NeXLCore.edgeenergy(z, ss, FFASTDB),
    ) : 0.0
ionizationcrosssection(z::Int, ss::Int, energy::AbstractFloat) =
    ionizationcrosssection(z, ss, energy, Bote2009)

"""
    jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) =

Compute the jump ratio.
"""
jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) = FFAST.jumpratio(z, ss)

include("strength.jl")

struct NeXL <: NeXLAlgorithm end

"""
    fluorescenceyield(z::Int, inner::Int, outer::Int)::Float64

The fraction of `inner` sub-shell ionizations that relax via a characteristic X-ray resulting from an
electronic transition from `outer` to `inner`.
"""
fluorescenceyield(z::Int, inner::Int, outer::Int, ::Type{NeXL})::Float64 =
    nexlTotalWeight(z, inner, inner, outer)

"""
    characteristicyield(z::Int, ionized::Int, inner::Int, outer::Int)::Float64

The fraction of `ionized` sub-shell ionizations that relax via a characteristic X-ray resulting from an
electronic transition from `outer` to `inner`.  This includes both direct transitions (where `outer`==`ionized`)
and cascade (where `outer` != `ionized` due to Coster-Kronig and previous decays.)
"""
characteristicyield(z::Int, ionized::Int, inner::Int, outer::Int, ::Type{NeXL})::Float64 =
    nexlTotalWeight(z, ionized, inner, outer)

"""
    characteristicXRayAvailable(z::Int, inner::Int, outer::Int)::Float64

Is the weight associated with this transition greater than zero?
"""
charactericXRayAvailable(z::Int, inner::Int, outer::Int, ::Type{NeXL})::Bool =
    nexlIsAvailable(z, inner, outer)
