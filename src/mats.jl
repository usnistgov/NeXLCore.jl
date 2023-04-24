"""
    Materials(
        name::AbstractString,
        els::AbstractArray{Element},
        ::Type{U},
        dims::Tuple;
        atomicweights::AbstractDict{Element,V} = Dict{Element,Float64}(),
        properties::AbstractDict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where {U<:AbstractFloat, V<:AbstractFloat}

A type to represent the composition of multiple points as a complement to `KRatios` and `HyperSpectrum`.
The data is stored in a much more memory efficient manner than `Array{Material}` would but can be accessed
either by coordinate like `mats[1,2]` returning a `Material` or by element like `mats[n"Fe"]` returning
an Array of mass fraction values.
"""
struct Materials{U<:AbstractFloat, V<:AbstractFloat, N} <: DenseArray{Material{U, V}, N}
    name::String
    planes::Dict{Element, Int}
    massfractions::Array{U}
    a::Dict{Element,V} # Optional: custom atomic weights for the keys in this Material
    properties::Dict{Symbol,Any} # :Density, :Description, :Pedigree, :Conductivity, ... + user defined

    function Materials(
        name::AbstractString,
        els::AbstractArray{Element},
        ::Type{U},
        dims::Tuple;
        atomicweights::AbstractDict{Element,V} = Dict{Element,Float64}(),
        properties::AbstractDict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where {U<:AbstractFloat, V<:AbstractFloat}
        plns = Dict( el => i for (i, el) in enumerate(els) )
        mfs = zeros(U, ( dims..., length(plns) ) )
        new{U, V, length(dims)}(name, plns, mfs, atomicweights, properties)
    end
end

function Base.setindex!(mats::Materials, mat::Material, idx::Int...)
    valof(mf::UncertainValue, ::Type{UncertainValue}) = mf
    valof(::Missing, ::Type{UncertainValue}) = zero(UncertainValue)
    valof(mf::UncertainValue, T::Type{<:AbstractFloat}) = T(value(mf))
    valof(mf::AbstractFloat, T::Type{<:AbstractFloat}) = T(mf)
    valof(::Missing, T::Type{<:AbstractFloat}) = zero(T)
    for (el, val) in mat.massfraction
        i = get(mats.planes, el, 0)
        if i > 0
            mats.massfractions[idx..., i] = valof(val, eltype(mats.massfractions))
        end
    end
end

Base.setindex!(mats::Materials, mat::Material, ci::CartesianIndex) = setindex!(mats, mat, ci.I...)

function Base.show(io::IO, mats::Materials)
    els = join(symbol.(sort( [ keys(mats.planes)... ])),", ")
    print(io, "Materials[$(mats.name), $(size(mats.massfractions)[1:end-1]), [$els]]")
end

function Base.getindex(mats::Materials{U,V,N}, elm::Element) where {U<:AbstractFloat,V<:AbstractFloat,N}
    i = get(mats.planes, elm, 0)
    i > 0 ? mats.massfractions[( 1:size(mats,i) for i in 1:N )..., i] : zeros(eltype(mat.massfractions), size(mat.massfractions)[1:end-1])
end

function Base.getindex(mats::Materials{U,V,N}, idx::Int...)::Material{U,V} where {U<:AbstractFloat,V<:AbstractFloat,N}
    mfs = Dict( el=>mats.massfractions[idx...,i] for (el, i) in mats.planes )
    Material("$(mats.name)[$idx]", mfs, mats.a, mats.properties)
end

Base.getindex(mats::Materials, ci::CartesianIndex) = getindex(mats, ci.I...)

function Base.getindex(mats::Materials{U,V,N}, idx::Int)::Material{U,V} where {U<:AbstractFloat,V<:AbstractFloat,N}
    mfs = Dict( el=>mats.massfractions[(idx-1)*size(mats.massfractions, ndims(mats.massfractions))+j] for (el, j) in mats.planes )
    Material("$(mats.name)[$idx]", mfs, mats.a, mats.properties)
end

Base.eachindex(mats::Materials) = CartesianIndices(size(mats))
Base.getindex(mats::Materials, sym::Symbol) = getindex(mats.properties, sym)
Base.setindex!(mats::Materials, val::Any, sym::Symbol) = setindex!(mats.properties, val, sym)
Base.CartesianIndices(mats::Materials)  = CartesianIndices(size(mats))
Base.size(mats::Materials) = size(mats.massfractions)[1:end-1]
Base.size(mats::Materials, dim::Int) = size(mats.massfractions, dim)
Base.ndims(mats::Materials) = ndims(mats.massfractions)-1
Base.eltype(::Materials{U,V,N}) where {U<:AbstractFloat,V<:AbstractFloat,N} = Material{U,V}
Base.similar(mats::Materials, T::Type=eltype(mats), dims=size(mats)) = Array{T}(undef, dims)

NeXLCore.a(elm::Element, mat::Materials) = get(mat.a, elm, a(elm))