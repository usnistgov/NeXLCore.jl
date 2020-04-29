"""
	minrequired(::Type{XXX})

Returns the minimum required properties.  Other classes implement this to check whether a Spectrum or Dict has all the
necessary properties for the specified algorithm or data structure.
"""
minproperties(::Type{Any}) =
    @assert "minproperties(...) needs to be specialized for this algorithm."

"""
    hasminrequired(ty::Type, item::Union{Spectrum,Dict{Symbol,Any}})

Does this spectrum have the minimal set of required properties?
"""
hasminrequired(ty::Type, dict::Dict{Symbol,Any}) = all(haskey(dict, a) for a in minproperties(ty))

"""
    requiredbutmissing(ty::Type, item::Union{Spectrum,Dict{Symbol,Any}})

List any required but missing properties.
"""
requiredbutmissing(ty::Type, dict::Dict{Symbol,Any}) = filter(a -> !haskey(dict, a), minproperties(ty))
