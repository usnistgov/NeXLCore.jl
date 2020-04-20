# Defines the Material struct and functions on it.
using Unitful
using Printf
using DataFrames

"""
    Material

Holds basic data about a material including name, density and composition in
mass fraction.

**Properties**

  - `:Density` Density in g/cm³
  - `:Description` Human friendly
  - `:Provenance` Source of data ("SRM-XXX", "Stoichiometry", "Wet-chemistry", "???")
  - `:Conductivity` "Insulator", "Semiconductor", "Conductor"
"""
struct Material
    name::String
    properties::Dict{Symbol,Any} # :Density, :Description, :Provenance
    a::Dict{Int,AbstractFloat} # Optional: custom atomic weights for the keys in this Material
    massfraction::Dict{Int,AbstractFloat}

"""
    Material(
        name::AbstractString,
        massfrac::Dict{Int,U},
        atomicweights::Dict{Int,V} = Dict{Int,Float64}(),
        properties::Dict{Symbol, Any} = Dict{Symbol, Any}()
    )
"""
    function Material(
        name::AbstractString,
        massfrac::Dict{Int,U},
        atomicweights::Dict{Int,V} = Dict{Int,Float64}(),
        properties::Dict{Symbol, Any} = Dict{Symbol, Any}()
    ) where {U<:AbstractFloat,V<:AbstractFloat}
        if sum(value.(values(massfrac)))>10.0
            @warn "The sum mass fraction is $(sum(values(massfrac))) which is much larger than unity."
        end
        new(name, properties, atomicweights, massfrac)
    end
end

elms(mat::Material) = Set(element(z) for z in keys(mat.massfraction))

import Base.*

function *(k::AbstractFloat, mat::Material)::Material
    mf = Dict{Int,AbstractFloat}((z, q * k) for (z, q) in mat.massfraction)
    return Material("$(k)×$(mat.name)", mf, mat.a, copy(mat.properties))
end

import Base.+

+(mat1::Material, mat2::Material)::Material = sum(mat1, mat2)

import Base.sum

function sum(
    mat1::Material,
    mat2::Material;
    name::Union{AbstractString, Missing}=missing,
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element,<:AbstractFloat}=Dict{Element,Float64}()
)::Material
    mf = Dict((elm, mat1[elm] + mat2[elm]) for elm in union(
        keys(mat1),
        keys(mat2),
    ))
    name = ismissing(name) ? "$(mat1.name)+$(mat2.name)" : name
    return material(name, mf, properties=properties, atomicweights=atomicweights)
end


"""
    name(mat::Material)

Return a human friendly short name for the Material.
"""
name(mat::Material) = mat.name

"""
    density(mat::Material)

Return the density in g/cm³ (Might be 'missing')
"""
density(mat::Material) = property(mat,:Density)
description(mat::Material) = property(mat,:Description)
provenance(mat::Material) = property(mat, :Provenance)

property(mat::Material, sym::Symbol) = get(mat.properties, sym, missing)


"""
    material(
        name::AbstractString,
        massfrac::Pair{Element,<:AbstractFloat}...;
        properties::Dict{Symbol,Any}=Dict{Symbol,Any)(),
        atomicweights::Dict{Element, <:AbstractFloat}=Dict{Element,Float64}()
    )
    material(
        name::AbstractString,
        massfrac::Dict{Element,<:AbstractFloat};
        properties::Dict{Symbol,Any}=Dict{Symbol,Any)(),
        atomicweights::Dict{Element, <:AbstractFloat}=Dict{Element,Float64}()
    )

Constuct a material from mass fraction pairs.
"""
function material(
    name::AbstractString,
    massfrac::Dict{Element,U};
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element, V}=Dict{Element,Float64}(),
) where {U <: AbstractFloat, V <: AbstractFloat}
    mf = Dict{Int,U}( (z(elm), v) for (elm, v) in massfrac)
    aw = Dict{Int,V}( (z(elm), v) for (elm, v) in atomicweights)
    return Material(name, mf, aw, properties)
end

material(
    name::AbstractString,
    massfrac::Pair{Element,<:AbstractFloat}...;
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element, <:AbstractFloat}=Dict{Element,Float64}(),
) = material(name, Dict(massfrac), atomicweights=atomicweights, properties=properties)


"""
    pure(elm::Element)

Construct a Material to represent a pure element.

Example:

    > pure(n"Fe")
"""
pure(elm::Element) =
    material("Pure $(symbol(elm))", Dict{}(elm=>1.0), properties=Dict{Symbol,Any}(:Density=>density(elm)))

function Base.show(io::IO, mat::Material)
    res="$(name(mat))["
    comma=""
    for (z, mf) in mat.massfraction
        res*=@sprintf("%s%s = %0.4f", comma, element(z).symbol, value(mf))
        comma=", "
    end
    if haskey(mat.properties,:Density)
        res*=@sprintf(", %0.2f g/cc)", density(mat))
    else
        res*="]"
    end
    print(io, res)
end

"""
    Base.convert(::Type{Material}, str::AbstractString)

Convert a DTSA-II style string into a material.
"""
function Base.convert(::Type{Material}, str::AbstractString)
    tmp = Dict{Element,Float64}()
    items = split(str,',')
    props = Dict{Symbol, Any}()
    name = item[1]
    for item in items[2:end]
        if (item[1]=='(') && (item[end]==')')
            zd = split(item[2:end-1],':')
            elm = parse(PeriodicTable.Element, zd[1])
            qty = parse(Float64,zd[2])/100.0
            tmp[elm]=qty
        else
            props[:Density] = parse(Float64,zd[2])
        end
    end
    return material(name, tmp, properties=props)
end


"""
    a(elm::Element, mat::Material)

Get the atomic weight for the specified Element in the specified Material.
"""
a(elm::Element, mat::Material) =
    get(mat.a, elm.number, a(elm))

Base.getindex(mat::Material, elm::Element) =
    get(mat.massfraction, elm.number, zero(eltype(values(mat.massfraction))))

Base.getindex(mat::Material, sym::Symbol) =
    property(mat, sym)

Base.get(mat::Material, sym::Symbol, def) =
    get(mat.properties, sym, def)

Base.setindex!(mat::Material, val, sym::Symbol) =
    mat.properties[sym] = val

nonneg(mat::Material, elm::Element) =
    max(0.0, value(mat[elm]))

"""
    normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}

Return the normalized mass fraction as a Dict{Element, AbstractFloat}.  Negative values
are set to zero.
"""
function normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}
    n = sum(nonneg(mat,elm) for elm in keys(mat))
    return Dict( (elm, nonneg(mat,elm)/n) for elm in keys(mat))
end

"""
    asnormalized(mat::Material, n=1.0)::Material

Convert the Material to a normalized Material form.  Negative mass fractions
are set to zero before normalization.
"""
function asnormalized(mat::Material, n=1.0)
    at = analyticaltotal(mat)
    if isapprox(at, n, rtol=1.0e-8) && startswith(name(mat),"N[")
        return mat
    else
        return Material(
            "N[$(name(mat)),$(n)]",
            Dict( (z(elm), n*nonneg(mat,elm)/at) for elm in keys(mat) ),
            mat.a,
            copy(mat.properties)
        )
    end
end

function Base.isapprox(mat1::Material, mat2::Material; atol=0.0, rtol=1.0e-4)
    for elm in union(keys(mat1),keys(mat2))
        if !isapprox(value(mat1[elm]), value(mat2[elm]), atol=atol, rtol=rtol)
            return false
        end
    end
    return true
end

"""
    massfraction(mat::Material)::Dict{Element, AbstractFloat}

The mass fraction as a Dict{Element, AbstractFloat}
"""
massfraction(mat::Material)::Dict{Element, AbstractFloat} =
    Dict( (element(z), mf) for (z, mf) in mat.massfraction )

"""
    keys(mat::Material)

Return an interator over the elements in the Material.
"""
Base.keys(mat::Material) =
    (element(z) for z in keys(mat.massfraction))

"""
    labeled(mat::Material)

Transform the mass fraction representation of a material into a Dict{MassFractionLabel,AbstractFloat}"""
labeled(mat::Material) =
    Dict( (MassFractionLabel(name(mat), element(z)), mf) for (z, mf) in mat.massfraction)

"""
    atomicfraction(mat::Material)::Dict{Element,AbstractFloat}

Return the composition in atomic fraction representation.
"""
function atomicfraction(mat::Material)::Dict{Element,AbstractFloat}
    norm = sum(mf/a(element(z),mat) for (z, mf) in mat.massfraction)
    return Dict( (element(z), (mf/a(element(z),mat))/norm) for (z, mf) in mat.massfraction)
end

"""
    analyticaltotal(mat::Material)

Return the sum of the positive mass fractions.
"""
analyticaltotal(mat::Material) =
    sum(nonneg(mat, elm) for elm in keys(mat))

"""
    haskey(mat::Material, elm::Element)

Does this material contain this element?
"""
Base.haskey(mat::Material, elm::Element) =
    haskey(mat.massfraction, z(elm))


"""
    atomicfraction(name::String, atomfracs::Pair{Element,Float64}...; properties::properties::Dict{Symbol, Any}, atomicweights::Dict{Element,Float64})

Build a Material from atomic fractions (or stoichiometries).
"""
function atomicfraction(
    name::AbstractString,
    atomfracs::Pair{Element,U}...;
    properties::Dict{Symbol, Any} = Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}()
)::Material where {U<:Real, V<:AbstractFloat}
    aw(elm) = get(atomicweights, elm, a(elm))
    norm = sum(af * aw(elm) for (elm, af) in atomfracs)
    massfracs = Dict((elm, (aw(elm) / norm) * af) for (elm, af) in atomfracs)
    return material(name, massfracs, atomicweights=atomicweights, properties=properties)
end

"""
    Base.parse(
        ::Type{Material},
        expr::AbstractString;
        name = missing,
        properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
        atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
    )::Material

Parse a Material from a string.
"""
function Base.parse(
    ::Type{Material},
    expr::AbstractString;
    name = missing,
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
)::Material where {V<:AbstractFloat}
    # Parses expressions like SiO2, Al2O3 or other simple (element qty) phrases
    function parseCompH1(expr::AbstractString)::Dict{PeriodicTable.Element, Int}
        parseSymbol(expr::AbstractString) =
            findfirst(z -> isequal(elements[z].symbol, expr), eachindex(elements))
        res = Dict{PeriodicTable.Element,Int}()
        start, idx = 1, collect(eachindex(expr))
        for i in eachindex(idx)
            if i<start
                continue
            elseif (i==start) || (i==start+1) # Abbreviations are 1 or 2 letters
                if (i == start) && !isuppercase(expr[idx[i]]) # Abbrevs start with cap
                    error("Element abbreviations must start with a capital letter. $(expr[idx[i]])")
                end
                next=i+1
                if (next>length(idx)) || isuppercase(expr[idx[next]]) || isdigit(expr[idx[next]])
                    z = parseSymbol(expr[idx[start]:idx[i]])
                    if isnothing(z)
                        error("Unrecognized element parsing compound: $(expr[start:i])")
                    end
                    elm, cx = elements[z], 1
                    if (next<=length(idx)) && isdigit(expr[idx[next]])
                        for stop in next:length(idx)
                            if (stop == length(idx)) || (!isdigit(expr[idx[stop+1]]))
                                cx = parse(Int, expr[idx[next]:idx[stop]])
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
                error("Can't interpret $(expr[idx[start]:idx[i]]) as an element.")
            end
        end
        return res
    end
    # Parses expressions like 'Ca5(PO4)3⋅(OH)'
    function parseCompH2(expr::AbstractString)::Dict{PeriodicTable.Element, Int}
        # println("Parsing: $(expr)")
        splt = split(expr, c->(c=='⋅') || (c=='.'))
        if length(splt)>1
            return mapreduce(parseCompH2, (a,b)->merge(+,a,b), splt)
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
            res, q = parseCompH2(expr[nextind(expr, start):prevind(expr,stop)]), 1
            if (nextind(expr,stop) > lastindex(expr)) || isdigit(expr[nextind(expr,stop)])
                for i in stop:lastindex(expr)
                    if (nextind(expr,i)>lastindex(expr)) || (!isdigit(expr[nextind(expr,i)]))
                        nxt = nextind(expr,stop)
                        q = nxt <= lastindex(expr) ? parse(Int, expr[nxt:i]) : 1
                        stop=i
                        break
                    end
                end
            end
            for elm in keys(res) res[elm] *= q end
            if start>1
                res=merge(+, res, parseCompH2(expr[1:prevind(expr,start)]))
            end
            if stop<lastindex(expr)
                res = merge(+, res,parseCompH2(expr[nextind(expr,stop):lastindex(expr)]))
            end
        else
            res = parseCompH1(expr)
        end
        # println("Parsing: $(expr) to $(res)")
        return res
    end

    # First split sums of Materials
    splt = split(expr, c -> c == '+')
    if length(splt) > 1
        return mapreduce(
            t -> parse(Material, strip(t)),
            (a, b) -> sum(a, b, name=name, properties=properties, atomicweights=atomicweights),
            splt,
        )
    end
    # Second handle N*Material where N is a number
    p = findfirst(c -> (c == '*') || (c == '×') || (c=='⋅'), expr)
    if !isnothing(p)
        return parse(Float64, expr[firstindex(expr):prevind(expr,p)]) *
               parse(Material, expr[nextind(expr,p):lastindex(expr)],
               name=expr[nextind(expr,p):lastindex(expr)], properties=properties, atomicweights=atomicweights)
    end
    # Then handle Material/N where N is a number
    p = findfirst(c -> c == '/', expr)
    if !isnothing(p)
        return (1.0 / parse(Float64, expr[nextidx(expr,p):lastindex(expr)])) *
               parse(Material, expr[firstindex(expr):prevind(expr,p)], name=expr[firstindex(expr):previdx(expr,p)], properties=properties, atomicweights=atomicweights)
    end
    # Finally parse material
    return atomicfraction(ismissing(name) ? expr : name, parseCompH2(expr)..., properties=properties, atomicweights=atomicweights)
end

macro mat_str(str)
    parse(Material,str)
end

"""
    NeXLUncertainties.asa(::Type{DataFrame}, mat::Material)

Tabulate the composition of this Material as a DataFrame.  Columns for
material name, element abbreviation, atomic number, atomic weight, mass fraction,
normalized mass fraction, and atomic fraction. Rows for each element in mat.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, mat::Material)
    res = DataFrame( Material = Vector{String}(), Element = Vector{String}(),
                AtomicNumber = Vector{Int}(), AtomicWeight = Vector{AbstractFloat}(),
                MassFraction = Vector{AbstractFloat}(), NormalizedMassFraction = Vector{AbstractFloat}(),
                AtomicFraction = Vector{AbstractFloat}() )
    af, tot = atomicfraction(mat), analyticaltotal(mat)
    for elm in sort(collect(keys(mat)))
        push!(res, ( name(mat), symbol(elm), z(elm), a(elm, mat), mat[elm], mat[elm]/tot, af[elm]) )
    end
    return res
end

"""
    NeXLUncertainties.asa(::Type{DataFrame}, mats::AbstractArray{Material}, mode=:MassFraction)

Tabulate the composition of a list of materials in a DataFrame.  One column
for each element in any of the materials.

    mode = :MassFraction | :NormalizedMassFraction | :AtomicFraction.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, mats::AbstractArray{Material}, mode=:MassFraction)
    elms = length(mats)==1 ? collect(keys(mats[1])) :
            Base.convert(Vector{Element}, sort(reduce(union, keys.(mats)))) # array of sorted Element
    cols = ( Symbol("Material"), Symbol.(symbol.(elms))..., Symbol("Total")) # Column names
    empty = NamedTuple{cols}( map(c->c==:Material ? Vector{String}() : Vector{AbstractFloat}(), cols) )
    res = DataFrame(empty) # Emtpy data frame with necessary columns
    for mat in mats
        vals = ( mode==:AtomicFraction ? atomicfraction(mat) :
                ( mode==:NormalizedMassFraction ? normalizedmassfraction(mat) :
                    massfraction(mat)))
        tmp = [ name(mat), (get(vals, elm, 0.0) for elm in elms)..., analyticaltotal(mat) ]
        push!(res, tmp)
    end
    return res
end

"""
    compare(unk::Material, known::Material)::DataFrame

Compare two compositions in a DataFrame.
"""
function compare(unk::Material, known::Material)::DataFrame
    um, km, z, kmf, rmf, dmf =
        Vector{String}(), Vector{String}(), Vector{String}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    fmf, kaf, raf, daf, faf =
        Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    afk, afr = atomicfraction(known), atomicfraction(unk)
    for elm in union(keys(known),keys(unk))
        push!(um, name(unk))
        push!(km, name(known))
        push!(z,elm.symbol)
        push!(kmf,known[elm])
        push!(rmf,unk[elm])
        push!(dmf,known[elm]-unk[elm])
        push!(fmf,(known[elm]-unk[elm])/known[elm])
        push!(kaf,get(afk, elm, 0.0))
        push!(raf,get(afr, elm, 0.0))
        push!(daf,get(afk, elm, 0.0)-get(afr, elm, 0.0))
        push!(faf,100.0*(get(afk, elm, 0.0)-get(afr, elm, 0.0))/get(afk, elm, 0.0))
    end
    return DataFrame(Unkown=um, Known=km, Elm=z, Cknown=kmf, Cresult=rmf, ΔC=dmf, ΔCoC=fmf, Aknown=kaf, Aresult=raf, ΔA=daf, ΔAoA=faf)
end

compare(unks::AbstractVector{Material}, known::Material) =
    mapreduce(unk->compare(unk, known), append!, unks)


"""
    mac(mat::Material, xray::Union{Float64,CharXRay}, alg::Type{<:NeXLAlgorithm}=FFASTDB)::Float64

Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, xray::Union{Float64,CharXRay}, alg::Type{<:NeXLAlgorithm}=FFASTDB) =
    mapreduce(elm->mac(elm, xray, alg)*mat[elm],+,keys(mat))


"""
    A structure defining a thin film of Material.
"""
struct Film
    material::Material
    thickness::AbstractFloat
end

Base.show(io::IO, flm::Film) =
    print(io, 1.0e7 * flm.thickness, " nm of ", name(flm.material))

"""
    transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat, alg::Type{<:NeXLAlgorithm}=FFASTDB) =

Compute the transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat, alg::Type{<:NeXLAlgorithm}=FFASTDB) =
    exp(-mac(flm.material, xrayE, alg) * csc(θ) * flm.thickness)

"""
    transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =

Compute the transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =
    transmission(flm, energy(cxr), θ)

material(film::Film) = film.material
thickness(film::Film) = film.thickness

function parsedtsa2comp(value::AbstractString)::Material
	try
		sp=split(value,",")
		name=sp[1]
		mf, props = Dict{Element,Float64}(), Dict{Symbol, Any}()
		for item in sp[2:end]
			if item[1]=='(' && item[end]==')'
				sp2=split(item[2:end-1],":")
				mf[parse(Element, sp2[1])]=0.01*parse(Float64,sp2[2])
			else
				props[:Density] = parse(Float64,item)
			end
		end
		return material(name, mf, properties=props)
	catch err
		warn("Error parsing composition $(value) - $(err)")
	end
end

function todtsa2comp(mat::Material)::String
    res=replace(name(mat),','=>'_')
    for (elm, qty) in mat.massfraction
        res*=",($(element(elm).symbol):$(100.0*qty))"
    end
    if haskey(mat.properties,:Density)
        res*=",$(mat.properties[:Density])"
    end
    return res
end

"""
    compositionlibrary()::Dict{String, Material}

Load the internal compositon library.
"""
function compositionlibrary()::Dict{String, Material}
    result = Dict{String, Material}()
    path = dirname(pathof(@__MODULE__))
    df = CSV.File("$(path)\\..\\data\\composition.csv") |> DataFrame
    for row in eachrow(df)
        name, density = row[1], row[2]
        elmc = collect(zip(element.(1:94), row[3:96])) # (i->getindex(row,i)).(3:96)))
        data=Dict{Element,Float64}(filter(a->(!ismissing(a[2]))&&(a[2]>0.0), elmc))
        properties=Dict{Symbol,Any}()
        if !ismissing(density)
            properties[:Density]=density
        end
        m = material(name, data, properties=properties)
        result[name] = m
    end
    return result
end
