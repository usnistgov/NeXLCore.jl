# Defines the Material struct and functions on it.
using Unitful
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
        density::Union{Missing,AbstractFloat},
        a::Dict{Int,V} = Dict{Int,Float64}(),
        description::Union{AbstractString,Missing} = missing,
    ) where {U<:AbstractFloat,V<:AbstractFloat}
        if sum(value.(values(massfrac)))>10.0
            @warn "The sum mass fraction is $(sum(values(massfrac))) which is much larger than unity."
        end
        props = Dict{Symbol,Any}()
        if (!ismissing(density)) && (density>0.0)
            props[:Density] = density
        end
        if !ismissing(description)
            props[:Description] = string(description)
        end
        new(name, props, a, massfrac)
    end
end

elms(mat::Material) = Set(element(z) for z in keys(mat.massfraction))

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

material(
    name::AbstractString,
    massfrac::Pair{Element,<:AbstractFloat}...;
    density::Union{Missing,AbstractFloat}=missing,
    a::Dict{Element, <:AbstractFloat}=Dict{Element,Float64}(),
    description::Union{Missing,AbstractString}=missing
) = material(name, Dict(massfrac), density, a, description)


"""
    pure(elm::Element)

Construct a Material to represent a pure element.

Example:

    > pure(n"Fe")
"""
pure(elm::Element) =
    material("Pure $(symbol(elm))", Dict{}(elm=>1.0), density(elm))

function Base.show(io::IO, mat::Material)
    res="$(name(mat))["
    comma=""
    for (z, mf) in mat.massfraction
        res*=@sprintf("%s%s = %0.4f", comma, element(z).symbol, value(mf))
        comma=", "
    end
    if !ismissing(density(mat))
        res*=@sprintf(", %0.2f g/cc)", density(mat))
    else
        res*="]"
    end
    print(io, res)
end

function Base.convert(::Type{Material}, str::AbstractString)
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

Base.get(mat::Material, sym::Symbol, def) =
    get(mat.properties, sym, def)


Base.setindex!(mat::Material, val, sym::Symbol) =
    mat.properties[sym] = val

nonneg(mat::Material, elm::Element) =
    max(0.0, value(mat[elm]))

"""
    normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}

The normalized mass fraction as a Dict{Element, AbstractFloat}.  Negative values
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
            mat[:Density],
            mat.a,
            mat[:Description]
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

The sum of the positive mass fractions.
"""
analyticaltotal(mat::Material) =
    sum(nonneg(mat,elm) for elm in keys(mat))

"""
    haskey(mat::Material, elm::Element)

Does this material contain this element?
"""
Base.haskey(mat::Material, elm::Element) =
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
    aw(elm) = get(atomicweights, elm, a(elm))
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
            (a, b) -> sum(a, b, name, density, atomicweights),
            splt,
        )
    end
    # Second handle N*Material where N is a number
    p = findfirst(c -> (c == '*') || (c == '×') || (c=='⋅'), expr)
    if !isnothing(p)
        return parse(Float64, expr[firstindex(expr):prevind(expr,p)]) *
               parse(Material, expr[nextind(expr,p):lastindex(expr)],
               name=expr[nextind(expr,p):lastindex(p)], density=density, atomicweights=atomicweights)
    end
    # Then handle Material/N where N is a number
    p = findfirst(c -> c == '/', expr)
    if !isnothing(p)
        return (1.0 / parse(Float64, expr[nextidx(expr,p):lastindex(expr)])) *
               parse(Material, expr[firstindex(expr):prevind(expr,p)], name=expr[firstindex(expr):previdx(expr,p)], density=density, atomicweights=atomicweights)
    end
    # Finally parse material
    return atomicfraction(ismissing(name) ? expr : name, parseCompH2(expr), density, atomicweights)
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

Compares two compositions in a DataFrame.

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
    mapreduce(elm->mac(elm, xray, alg)*massfraction(elm,mat),+,keys(mat))


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

Transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, xrayE::AbstractFloat, θ::AbstractFloat, alg::Type{<:NeXLAlgorithm}=FFASTDB) =
    exp(-mac(flm.material, xrayE, alg) * csc(θ) * flm.thickness)

"""
    transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =

Transmission fraction of an X-ray at the specified angle through a Film.
"""
transmission(flm::Film, cxr::CharXRay, θ::AbstractFloat) =
    transmission(flm, energy(cxr), θ)

material(film::Film) = film.material
thickness(film::Film) = film.thickness



function parsedtsa2comp(value::AbstractString)::Material
	try
		sp=split(value,",")
		name=sp[1]
		mf = Dict{Element,Float64}()
		den = missing
		for item in sp[2:end]
			if item[1]=='(' && item[end]==')'
				sp2=split(item[2:end-1],":")
				mf[parse(Element, sp2[1])]=0.01*parse(Float64,sp2[2])
			else
				den = parse(Float64,item)
			end
		end
		return material(name, mf, den)
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

function compositionlibrary()::Dict{String, Material}
    result = Dict{String, Material}()
    path = dirname(pathof(@__MODULE__))
    df = CSV.File("$(path)\\..\\data\\composition.csv") |> DataFrame
    for row in eachrow(df)
        name, density = row[1], row[2]
        elmc = collect(zip(element.(1:94), row[3:96])) # (i->getindex(row,i)).(3:96)))
        data=Dict{Element,Float64}(filter(a->(!ismissing(a[2]))&&(a[2]>0.0), elmc))
        m = material(name, data, density)
        result[name] = m
    end
    return result
end
