import BoteSalvatICX # For ionization crosssections

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
        edgeenergy(z, ss),
    ) : 0.0