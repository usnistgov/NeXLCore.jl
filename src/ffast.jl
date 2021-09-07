import FFAST # for mass absorption coefficienct

"""
FFAST represents an implementation of mass-absorption coefficients and associated
edge and other elemental data.  It is a thin wrapper around the FFAST library.
"""
struct FFASTDB <: NeXLAlgorithm end

mac(elm::Element, energy::Float64, ::Type{FFASTDB})::Float64 =
    FFAST.mac(FFAST.PhotoElectricMAC, z(elm), energy)

function macU(elm::Element, energy::Float64, ::Type{FFASTDB})::UncertainValue
    macv = FFAST.mac(FFAST.PhotoElectricMAC, z(elm), energy)
    return uv(
        macv,
        min(FFAST.fractionaluncertainty(FFAST.SolidLiquid, z(elm), energy)[1], 0.9) * macv,
    )
end

edgeenergy(z::Int, ss::Int, ::Type{FFASTDB})::Float64 = FFAST.edgeenergy(z, ss)

eachelement(::Type{FFASTDB}) = FFAST.eachelement()

subshellsindexes(z::Int, ::Type{FFASTDB}) = FFAST.eachedge(z)

"""
    energy(z::Int, inner::Int, outer::Int)::Float64

Return energy (in eV) of the transition by specified inner and outer sub-shell index.
"""
energy(z::Int, inner::Int, outer::Int, ::Type{FFASTDB})::Float64 =
    FFAST.edgeenergy(z, inner) - FFAST.edgeenergy(z, outer)

"""
    jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) =

Compute the jump ratio.
"""
jumpratio(z::Int, ss::Int, ::Type{FFASTDB}) = FFAST.jumpratio(z, ss)
