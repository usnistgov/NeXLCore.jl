# Defines the Material struct and functions on it.
using Unitful
using PeriodicTable
using Printf
using DataFrames

"""
    Material

Holds basic data about a material including name, density and composition in
mass fraction.
"""
struct Material
    name::String
    properties::Dict{Symbol,Any} # :Density, :Description,
    a::Dict{Int,AbstractFloat} # Optional: custom atomic weights for the keys in this Material
    massfraction::Dict{Int,AbstractFloat}
    function Material(
        name::AbstractString,
        massfrac::Dict{Int,U},
        density::AbstractFloat,
        a::Dict{Int,V} = Dict{Int,Float64}(),
        description::Union{AbstractString,Missing} = missing,
    ) where {U<:AbstractFloat,V<:AbstractFloat}
        props = Dict{Symbol,Any}()
        if !ismissing(density)
            props[:Density] = density
        end
        if !ismissing(description)
            props[:Description] = string(description)
        end
        new(name, props, a, massfrac)
    end
end

import Base.*

function *(k::AbstractFloat, mat::Material)::Material
    mf = Dict{Int,AbstractFloat}((z, q * k) for (z, q) in mat.massfraction)
    return Material("$(k)×$(mat.name)", mf, density(mat), mat.a)
end

import Base.+

+(mat1::Material, mat2::Material)::Material = sum(mat1, mat2, missing, missing)

import Base.sum

function sum(
    mat1::Material,
    mat2::Material,
    name::Union{Missing,AbstractString} = missing,
    density::Union{Missing,AbstractFloat} = missing,
    atomicweight::Dict{Element,<:AbstractFloat} = Dict{Element,Float64}()
)::Material
    mf = Dict{Element,AbstractFloat}((elm, mat1[elm] + mat2[elm]) for elm in union(
        keys(mat1),
        keys(mat2),
    ))
    return material(
        ismissing(name) ? "$(mat1.name) + $(mat2.name)" : name,
        mf,
        density,
        atomicweight
    )
end


"""
    name(mat::Material)

Human friendly short name for the Material.
"""
name(mat::Material) = mat.name

"""
    density(mat::Material)

Density in g/cm³ (Might be 'missing')
"""
density(mat::Material) = property(mat,:Density)

description(mat::Material) = property(mat,:Description)

property(mat::Material, sym::Symbol) = get(mat.properties, sym, missing)

"""
    material(
      name::AbstractString,
      massfrac::Dict{Element,<:AbstractFloat},
      density=missing,
      a::Dict{Element,<:AbstractFloat}=Dict{Element,Float64}(),
      description::Union{Missing,AbstractString}=missing
    )

Construct a material from mass fractions and (optional) atomic weights.
"""
function material(name::AbstractString,
    massfrac::Dict{Element, U},
    density::Union{Missing,AbstractFloat}=missing,
    a::Dict{Element, V}=Dict{Element,Float64}(),
    description::Union{Missing,AbstractString}=missing
) where { U <: AbstractFloat, V <: AbstractFloat }
    mf = Dict( (z(elm), val) for (elm, val) in massfrac )
    aw = Dict( (z(elm), val) for (elm, val) in a )
    return Material(name, mf, density, aw, description)
end

"""
    pure(elm::Element)

Construct a Material to represent a pure element.

Example:

    > pure(n"Fe")
"""
pure(elm::Element) =
    material("Pure "*symbol(elm), Dict{}(elm=>1.0), density(elm))

function Base.show(io::IO, mat::Material)
    print(io, name(mat)," = (")
    comma=""
    for (z, mf) in mat.massfraction
        print(io, @sprintf("%s%s = %0.4f", comma, element(z).symbol, mf))
        comma=", "
    end
    if !ismissing(density(mat))
        print(io, @sprintf(", %0.2f g/cc", density(mat)))
    end
    print(io,")")
end

function convert(::Type{Material}, str::AbstractString)
    tmp = Dict{Element,Float64}()
    items = split(str,',')
    density = missing
    name = item[1]
    for item in items[2:end]
        if (item[1]=='(') && (item[end]==')')
            zd = split(item[2:end-1],':')
            elm = parse(PeriodicTable.Element,zd[1])
            qty = parse(Float64,zd[2])/100.0
            tmp[elm]=qty
        else
            density = parse(Float64,zd[2])
        end
    end
    return Material(name, tmp, density)
end


"""
    a(elm::Element, mat::Material)

Get the atomic weight for the specified Element in the specified Material.
"""
a(elm::Element, mat::Material) =
    get(mat.a, elm.number, a(elm))

massfraction(elm::Element, mat::Material) =
    get(mat.massfraction,elm.number,zero(eltype(values(mat.massfraction))))

Base.getindex(mat::Material, elm::Element) =
    massfraction(elm,mat)

Base.getindex(mat::Material, sym::Symbol) =
    property(mat, sym)

"""
    normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}

The normalized mass fraction as a Dict{Element, AbstractFloat}
"""
function normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}
    n = sum(values(mat.massfraction))
    return Dict( (element(z), mf/n) for (z, mf) in mat.massfraction )
end

"""
    massfraction(mat::Material)::Dict{Element, AbstractFloat}

The mass fraction as a Dict{Element, AbstractFloat}
"""
massfraction(mat::Material)::Dict{Element, AbstractFloat} =
    Dict( (element(z), mf) for (z, mf) in mat.massfraction )

"""
    keys(mat::Material)

Returns an interator over the elements in the Material.
"""
Base.keys(mat::Material) =
    (element(z) for z in keys(mat.massfraction))

"""
    labeled(mat::Material)

Transform the mass fraction representation of a material into a Dict{MassFractionLabel,AbstractFloat}"""
labeled(mat::Material) =
    Dict( (MassFractionLabel(element(z), mat), mf) for (z, mf) in mat.massfraction)

"""
    atomicfraction(mat::Material)::Dict{Element,AbstractFloat}

The composition in atomic fraction representation.
"""
function atomicfraction(mat::Material)::Dict{Element,AbstractFloat}
    norm = sum(mf/a(element(z),mat) for (z, mf) in mat.massfraction)
    return Dict( (element(z), (mf/a(element(z),mat))/norm) for (z, mf) in mat.massfraction)
end

"""
    analyticaltotal(mat::Material)

The sum of the mass fractions.
"""
analyticaltotal(mat::Material) =
    sum(values(mat.massfraction))

"""
    has(mat::Material, elm::Element)

Does this material contain this element?
"""
has(mat::Material, elm::Element) =
    haskey(mat.massfraction, z(elm))


"""
    atomicfraction(name::String, atomfracs::Dict{Element,Float64}, density = nothing, atomicweights::Dict{Element,Float64}=Dict())

Build a Material from atomic fractions (or stoichiometries).
"""
function atomicfraction(
    name::AbstractString,
    atomfracs::Dict{Element,U},
    density::Union{Missing,AbstractFloat} = missing,
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
)::Material where {U<:Number,V<:AbstractFloat}
    aw(elm) = get(atomicweights, elm.number, a(elm))
    norm = sum(af * aw(elm) for (elm, af) in atomfracs)
    massfracs = Dict((elm, af * aw(elm) / norm) for (elm, af) in atomfracs)
    return material(name, massfracs, density, atomicweights)
end

function Base.parse(
    ::Type{Material},
    expr::AbstractString;
    name = missing,
    density = missing,
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
)::Material where {V<:AbstractFloat}
    # First split sums of Materials
    tmp = split(expr, c -> c == '+')
    if length(tmp) > 1
        return mapreduce(
            t -> parse(Material, strip(t)),
            (a, b) -> sum(a, b, name, density, atomicweights),
            tmp,
        )
    end
    # Second handle N*Material where N is a number
    p = findfirst(c -> (c == '*') || (c == '×'), expr)
    if !isnothing(p)
        return parse(Float64, expr[1:p-1]) *
               parse(Material, expr[p+1:end],
               name=expr[p+1:end], density=density,
               atomicweights=atomicweights)
    end
    # Then handle Material/N where N is a number
    p = findfirst(c -> c == '/', expr)
    if !isnothing(p)
        return (1.0 / parse(Float64, expr[p+1:end])) *
               parse(Material, expr[1:p-1], name=expr[1:p-1], density=density, atomicweights=atomicweights)
    end
    # Finally parse material
    return atomicfraction(ismissing(name) ? expr : name, parseCompH2(expr), density, atomicweights)
end

# Parses expressions like 'Ca5(PO4)3⋅(OH)'
function parseCompH2(expr::AbstractString)::Dict{PeriodicTable.Element, Int}
    # println("Parsing: $(expr)")
    tmp = split(expr, c->(c=='⋅') || (c=='.'))
    if length(tmp)>1
        println(tmp)
        return mapreduce(parseCompH2, merge, tmp)
    end
    cx, start, stop = 0, -1, -1
    for i in eachindex(expr)
        if expr[i]=='('
            if cx==0 start=i end
            cx+=1
        end
        if expr[i]==')'
            cx-=1
            if cx<0
                error("Unmatched right parenthesis.")
            elseif cx==0
                stop = i
                break
            end
        end
    end
    if (start>0) && (stop>start)
        tmp, q = parseCompH2(expr[start+1:stop-1]), 1
        if (stop+1 > length(expr)) || isdigit(expr[stop+1])
            for i in stop:length(expr)
                if (i+1>length(expr)) || (!isdigit(expr[i+1]))
                    q = parse(Int, expr[stop+1:i])
                    stop=i
                    break
                end
            end
        end
        for elm in keys(tmp) tmp[elm] *= q end
        if start>1
            merge!(tmp,parseCompH2(expr[1:start-1]))
        end
        if stop<length(expr)
            merge!(tmp,parseCompH2(expr[stop+1:end]))
        end
    else
        tmp = parseCompH1(expr)
    end
    # println("Parsing: $(expr) to $(tmp)")
    return tmp
end


# Parses expressions like SiO2, Al2O3 or other simple (element qty) phrases
function parseCompH1(expr::AbstractString)::Dict{PeriodicTable.Element,Int}
    parseSymbol(expr::AbstractString) =
        findfirst(z -> isequal(elements[z].symbol, expr), 1:length(elements))
    res = Dict{PeriodicTable.Element,Int}()
    start = 1
    for i in eachindex(expr)
        if i<start
            continue
        elseif (i==start) || (i==start+1) # Abbreviations are 1 or 2 letters
            if (i == start) && !isuppercase(expr[i]) # Abbrevs start with cap
                error("Element abbreviations must start with a capital letter. $(expr[i])")
            end
            next=i+1
            if (next>length(expr)) || isuppercase(expr[next]) || isdigit(expr[next])
                z = parseSymbol(expr[start:i])
                if isnothing(z)
                    error("Unrecognized element parsing compound: $(expr[start:i])")
                end
                elm, cx = elements[z], 1
                if (next<=length(expr)) && isdigit(expr[next])
                    for stop in next:length(expr)
                        if (stop == length(expr)) || (!isdigit(expr[stop+1]))
                            cx = parse(Int, expr[next:stop])
                            start = stop + 1
                            break
                        end
                    end
                else
                    start = next
                end
                res[elm] = get(res, elm, 0) + cx
            end
        else
            error("Can't interpret $(expr[start:i]) as an element.")
        end
    end
    return res
end


"""
    summarize(mat::Material)

Summarize the composition of this Material as a DataFrame.  Columns for
material name, element abbreviation, atomic number, atomic weight, mass fraction,
normalized mass fraction, and atomic fraction. Rows for each element in mat.
"""
function summarize(mat::Material)
    res = DataFrame( Material = Vector{String}(), Element = Vector{String}(),
                AtomicNumber = Vector{Int}(), AtomicWeight = Vector{AbstractFloat}(),
                MassFraction = Vector{AbstractFloat}(), NormalizedMassFraction = Vector{AbstractFloat}(),
                AtomicFraction = Vector{AbstractFloat}() )
    af, tot = atomicfraction(mat), analyticaltotal(mat)
    for elm in sort(collect(keys(mat)))
        push!(res, ( name(mat), symbol(elm), z(elm), a(elm), mat[elm], mat[elm]/tot, af[elm]) )
    end
    return res
end


"""
    summarize(mats::AbstractArray{Material}, mode=:MassFraction)

Summarize the composition of a list of materials in a DataFrame.  One column
for each element in any of the materials.

    mode = :MassFraction | :NormalizedMassFraction | :AtomicFraction.
"""
function summarize(mats::AbstractArray{Material}, mode=:MassFraction)
    elms = convert(Array{Element}, sort(reduce(union, keys.(mats)))) # array of sorted Element
    cols = ( Symbol("Material"), Symbol.(symbol.(elms))...) # Column names
    empty = NamedTuple{cols}( map(c->c==:Material ? Vector{String}() : Vector{AbstractFloat}(), cols) )
    res = DataFrame(empty) # Emtpy data frame with necessary columns
    for mat in mats
        vals = ( mode==:AtomicFraction ? atomicfraction(mat) :
                ( mode==:NormalizedMassFraction ? normalizedmassfraction(mat) :
                    massfraction(mat)))
        tmp = [ name(mat), (get(vals, elm, 0.0) for elm in elms)... ]
        push!(res, tmp)
    end
    return res
end


"""
    mac(mat::Material, xray::Union{Float64,CharXRay})::Float64

Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, xray::Union{Float64,CharXRay}) =
    mapreduce(elm->mac(elm, xray)*massfraction(elm,mat),+,keys(mat))


"""
    A structure defining a thin film of Material.
"""
struct Film
    material::Material
    thickness::AbstractFloat
end

Base.show(io::IO, flm::Film) =
    print(io, 1.0e7 * flm.thickness, " nm of ", name(flm.coating))

"""
    transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat) =

Transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat) =
    exp(-mac(flm.material, xrayE) * csc(θ) * flm.thickness)

"""
    transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =

Transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =
    transmission(flm, energy(cxr), θ)
