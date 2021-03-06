# Defines the Material struct and functions on it.
using Unitful
using Printf
using DataFrames
using LaTeXStrings
import Base.rand
import Statistics

"""
    Material

Holds basic data about a material including name, composition in mass fraction and optional propreties.

**Properties**

  - `:Density` Density in g/cm³
  - `:Description` Human friendly
  - `:Pedigree` Quality indicator for compositional data ("SRM-XXX", "CRM-XXX", "NIST K-Glass", "Stoichiometry", "Wet-chemistry by ???", "WDS by ???", "???")
  - `:Conductivity` => :Insulator, :Semiconductor, :Conductor
"""
struct Material{U<:AbstractFloat,V<:AbstractFloat}
    name::String
    properties::Dict{Symbol,Any} # :Density, :Description, :Pedigree, :Conductivity, ... + user defined
    massfraction::Dict{Int,U}
    a::Dict{Int,V} # Optional: custom atomic weights for the keys in this Material
    
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
        properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where {U<:AbstractFloat, V<:AbstractFloat}
        if sum(value.(values(massfrac))) > 10.0
            @warn "The sum mass fraction is $(sum(values(massfrac))) which is much larger than unity."
        end
        new{U,V}(name, properties, massfrac, atomicweights)
    end
end

rename(mat::Material, newname::AbstractString) = Material(newname, mat.massfraction, mat.a, mat.properties)

Base.copy(m::Material) =
    Material(m.name, copy(m.massfraction), copy(m.a), copy(m.properties))

elms(mat::Material) = Set(element(z) for z in keys(mat.massfraction))

function Base.:*(k::AbstractFloat, mat::Material)::Material
    mf = Dict((z, q * k) for (z, q) in mat.massfraction)
    return Material("$(k)×$(mat.name)", mf, copy(mat.a), copy(mat.properties))
end

Base.isequal(m1::Material, m2::Material) =
    isequal(m1.name, m2.name) &&
    isequal(m1.properties, m2.properties) &&
    isequal(m1.a, m2.a) &&
    isequal(m1.massfraction, m2.massfraction)

Base.hash(m::Material, h::UInt) =
    hash(m.name, hash(m.properties, hash(m.a, hash(m.massfraction, h))))


Base.:+(mat1::Material, mat2::Material)::Material = sum(mat1, mat2)


"""
    Base.sum(
        mat1::Material,
        mat2::Material;
        name::Union{AbstractString,Missing} = missing,
        properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
        density::Union{Missing,AbstractFloat} = missing,
        description::Union{Missing,AbstractString} = missing,
        pedigree::Union{Missing,AbstractString} = missing,
        conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
    )::Material


Construct a Material that represents the mass-fraction sum of mat1 and mat2. This function 
is often used along with Base.:*(k::AbstractFloat, mat::Material)::Material to construct
mixtures of compounds.  Ultimately, expressions like `mat"0.5*Al2O3+0.5*MgO"` are computed
using `sum(...)`.
"""
function Base.sum(
    mat1::Material,
    mat2::Material;
    name::Union{AbstractString,Missing} = missing,
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
)::Material
    plus(v1::AbstractFloat, v2::AbstractFloat) = 
        (σ(v1)==0.0)&&(σ(v2)==0.0) ? value(v1)+value(v2) : uv(value(v1)+value(v2),sqrt(σ(v1)^2+σ(v2)^2))
    mf = Dict((elm, plus(mat1[elm], mat2[elm])) for elm in union(keys(mat1), keys(mat2)))
    aw = Dict{Element,Float64}()
    for z in union(keys(mat1.a), keys(mat2.a))
        elm = elements[z]
        aw[elm] =
            (value(mat1[elm]) + value(mat2[elm])) / (value(mat1[elm]) / a(elm, mat1) + value(mat2[elm]) / a(elm, mat2))
    end
    name = ismissing(name) ? "$(mat1.name)+$(mat2.name)" : name
    return material(
        name,
        mf;
        properties = properties,
        atomicweights = aw,
        density = density,
        description = description,
        conductivity = conductivity,
        pedigree = pedigree,
    )
end

"""
    Base.sum(
      data::Dict{Material, <:AbstractFloat};
      name::Union{AbstractString,Missing} = missing,
      properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
      density::Union{Missing,AbstractFloat} = missing,
      description::Union{Missing,AbstractString} = missing,
      pedigree::Union{Missing,AbstractString} = missing,
      conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
)::Material

Sum together proportions of various `Material` structs.  The dictionary defines the material and the mass fraction of that material.
"""
function Base.sum(
    data::Dict{Material, <:AbstractFloat};
    name::Union{AbstractString,Missing} = missing,
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
)::Material
    assign(val, prop, props) = if !ismissing(val) props[prop]=val end
    res = mapreduce((a,b)->sum(a,b), data) do (mat, f)
        f*mat
    end
    props = copy(properties)
    assign(density, :Density, props)
    assign(description, :Description, props)
    assign(pedigree, :Pedigree, props)
    assign(conductivity, :Conductivity, props)
    return Material(ismissing(name) ? res.name : name, res.massfraction, res.a, props)
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
density(mat::Material) = property(mat, :Density)
description(mat::Material) = property(mat, :Description)
pedigree(mat::Material) = property(mat, :Pedigree)

"""
    atoms_per_cm³(mat::Material, elm::Element) =

Number of atoms per cm³ of the specified Element in the specified Material.  The Material must define
the `:Density` property.

    atoms_per_cm³(mat::Material)

Total number of atoms per cm³ for all elements in mat. 
"""
atoms_per_cm³(mat::Material, elm::Element) =
    mat[:Density] * mat[elm] * 6.0221366516752e23 / a(elm, mat)
atoms_per_cm³(mat::Material) = sum(atoms_per_cm³(mat, elm) for elm in keys(mat))

property(mat::Material, sym::Symbol) = get(mat.properties, sym, missing)


"""
    material(
        name::AbstractString,
        massfrac::Pair{Element,<:AbstractFloat}...;
        properties::Dict{Symbol,Any}=Dict{Symbol,Any)(),
        atomicweights::Dict{Element, <:AbstractFloat}=Dict{Element,Float64}(),
        density::Union{Missing, AbstractFloat}=missing,
        description::Union{Missing, AbstractString}=missing,
        pedigree::Union{Missing, AbstractString}=missing,
        conductivity::Union{Missing, Symbol}=missing, # :Conductor, :Semiconductor, :Insulator
    )
    material(
        name::AbstractString,
        massfrac::Dict{Element,<:AbstractFloat};
        properties::Dict{Symbol,Any}=Dict{Symbol,Any)(),
        atomicweights::Dict{Element, <:AbstractFloat}=Dict{Element,Float64}(),
        density::Union{Missing, AbstractFloat}=missing,
        description::Union{Missing, AbstractString}=missing,
        pedigree::Union{Missing, AbstractString}=missing,
        conductivity::Union{Missing, Symbol}=missing, # :Conductor, :Semiconductor, :Insulator
    )

Constuct a material from mass fraction pairs.
"""
function material(
    name::AbstractString,
    massfrac::Dict{Element,U};
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
) where {U<:AbstractFloat,V<:AbstractFloat}
    props = copy(properties)
    (!ismissing(density)) && ((props[:Density] = density) == density)
    (!ismissing(description)) && ((props[:Description] = description) == description)
    (!ismissing(pedigree)) && ((props[:Pedigree] = pedigree) == pedigree)
    (!ismissing(conductivity)) && ((props[:Conductivity] = conductivity) == conductivity)
    reallyhas(kwa, sym) =
        haskey(kwa, sym) && (!ismissing(kwa[sym])) && (!isnothing(kwa[sym]))
    mf = Dict{Int,U}((z(elm), v) for (elm, v) in massfrac)
    aw = Dict{Int,V}((z(elm), v) for (elm, v) in atomicweights)
    return Material(name, mf, aw, props)
end

material(
    name::AbstractString,
    massfrac::Pair{Element,U}...;
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
) where {U<:AbstractFloat,V<:AbstractFloat} = material(
    name,
    Dict(massfrac);
    properties = properties,
    atomicweights = atomicweights,
    density = density,
    description = description,
    pedigree = pedigree,
    conductivity = conductivity,
)


"""
     material(str::String, density::Float64)

Similar to `mat"..."` except requires you to specify a density.
"""
material(str::String, density::Float64) = parse(Material, str, density = density)

"""
    pure(elm::Element)

Construct a Material to represent a pure element.

Example:

    > pure(n"Fe")
"""
pure(elm::Element) =
    material("Pure $(symbol(elm))", Dict{}(elm => 1.0), density = density(elm))

function Base.show(io::IO, mat::Material)
    res = "$(name(mat))["
    res *= join(
        (
            @sprintf("%s=%0.4f", element(z).symbol, value(mf)) for
            (z, mf) in mat.massfraction
        ),
        ",",
    )
    if haskey(mat.properties, :Density)
        res *= @sprintf(",%0.2f g/cm³", density(mat))
    end
    res *= "]"
    print(io, res)
end

"""
    Base.convert(::Type{Material}, str::AbstractString)

Convert a DTSA-II style string into a material.
"""
function Base.convert(::Type{Material}, str::AbstractString)
    tmp = Dict{Element,Float64}()
    items = split(str, ',')
    density = missing
    name = item[1]
    for item in items[2:end]
        if (item[1] == '(') && (item[end] == ')')
            zd = split(item[2:end-1], ':')
            elm = parse(PeriodicTable.Element, zd[1])
            qty = parse(Float64, zd[2]) / 100.0
            tmp[elm] = qty
        else
            density = parse(Float64, zd[2])
        end
    end
    return material(name, tmp, density = density)
end


"""
    a(elm::Element, mat::Material)

Get the atomic weight for the specified Element in the specified Material.
"""
a(elm::Element, mat::Material) = get(mat.a, elm.number, a(elm))

Base.getindex(mat::Material, elm::Element) =
    get(mat.massfraction, elm.number, zero(eltype(values(mat.massfraction))))
Base.getindex(mat::Material, z::Int) =
    get(mat.massfraction, z, zero(eltype(values(mat.massfraction))))

Base.getindex(mat::Material, sym::Symbol) = property(mat, sym)

Base.get(mat::Material, sym::Symbol, def) = get(mat.properties, sym, def)

Base.setindex!(mat::Material, val, sym::Symbol) = mat.properties[sym] = val

nonneg(mat::Material, elm::Element) = max(0.0, value(mat[elm]))

"""
    normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}

Return the normalized mass fraction as a Dict{Element, AbstractFloat}.  Negative values
are set to zero.
"""
function normalizedmassfraction(mat::Material)::Dict{Element,AbstractFloat}
    n = analyticaltotal(mat)
    return Dict((elm, nonneg(mat, elm) / n) for elm in keys(mat))
end

normalized(mat::Material, elm::Element) = nonneg(mat, elm) / analyticaltotal(mat)
"""
    asnormalized(mat::Material, n=1.0)::Material

Convert the Material to a normalized Material form.  Negative mass fractions
are set to zero before normalization.
"""
function asnormalized(mat::Material, n = 1.0)
    at = analyticaltotal(mat)
    if isapprox(at, n, rtol = 1.0e-8) && startswith(name(mat), "N[")
        return mat
    else
        return Material(
            "N[$(name(mat)),$(n)]",
            Dict((z(elm), n * nonneg(mat, elm) / at) for elm in keys(mat)),
            mat.a,
            copy(mat.properties),
        )
    end
end


"""
    Base.isapprox(mat1::Material, mat2::Material; atol = 1.0e-4)

Are these Material(s) equivalent to within `atol`?
"""
function Base.isapprox(mat1::Material, mat2::Material; atol = 1.0e-4)
    return all(
        isapprox(value(mat1[elm]), value(mat2[elm]), atol = atol) && #
        isapprox(σ(mat1[elm]), σ(mat2[elm]), atol = atol) 
            for elm in union(keys(mat1), keys(mat2))
    )
end

"""
    massfraction(mat::Material)::Dict{Element, AbstractFloat}

The mass fraction as a Dict{Element, AbstractFloat}
"""
massfraction(mat::Material)::Dict{Element,AbstractFloat} =
    Dict((element(z), mf) for (z, mf) in mat.massfraction)

"""
    keys(mat::Material)

Return an interator over the elements in the Material.
"""
Base.keys(mat::Material) = (element(z) for z in keys(mat.massfraction))

"""
    labeled(mat::Material)

Transform the mass fraction representation of a material into a Dict{MassFractionLabel,AbstractFloat}"""
labeled(mat::Material) =
    Dict((MassFractionLabel(name(mat), element(z)), mf) for (z, mf) in mat.massfraction)

"""
    atomicfraction(mat::Material)::Dict{Element,AbstractFloat}

Return the composition in atomic fraction representation.
"""
function atomicfraction(mat::Material)::Dict{Element,AbstractFloat}
    norm = sum(mf / a(element(z), mat) for (z, mf) in mat.massfraction)
    return Dict(
        (element(z), (mf / a(element(z), mat)) / norm) for (z, mf) in mat.massfraction
    )
end

"""
    analyticaltotal(mat::Material)

Return the sum of the positive mass fractions.
"""
analyticaltotal(mat::Material) = sum(nonneg(mat, elm) for elm in keys(mat))

"""
    haskey(mat::Material, elm::Element)

Does this material contain this element?
"""
Base.haskey(mat::Material, elm::Element) = haskey(mat.massfraction, z(elm))

"""
    haskey(mat::Material, sym::Symbol)

Does this material have this property defined?
"""
Base.haskey(mat::Material, sym::Symbol) = haskey(mat.properties, sym)


"""
    atomicfraction(
        name::String,
        atomfracs::Union{Dict{Element,Float64},Pair{Element,Float64}...};
        properties::properties::Dict{Symbol, Any},
        atomicweights::Dict{Element,Float64},
        density::Union{Missing, AbstractFloat}=missing,
        description::Union{Missing, AbstractString}=missing,
        pedigree::Union{Missing, AbstractString}=missing,
        conductivity::Union{Missing, Symbol}=missing, # :Conductor, :Semiconductor, :Insulator
) # density, description, pedigree, conductivity

Build a Material from atomic fractions (or stoichiometries).
"""
atomicfraction(
    name::AbstractString,
    atomfracs::Pair{Element,U}...;
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
) where {U<:Real,V<:AbstractFloat} = atomicfraction(
    name,
    Dict(atomfracs);
    properties = properties,
    atomicweights = atomicweights,
    density = density,
    description = description,
    pedigree = pedigree,
    conductivity = conductivity,
)

function atomicfraction(
    name::AbstractString,
    atomfracs::Dict{Element,U};
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
) where {U<:Real,V<:AbstractFloat}
    aw(elm) = get(atomicweights, elm, a(elm))
    norm = sum(af * aw(elm) for (elm, af) in atomfracs)
    massfracs = Dict((elm, (aw(elm) / norm) * af) for (elm, af) in atomfracs)
    return material(
        name,
        massfracs;
        atomicweights = atomicweights,
        properties = properties,
        density = density,
        description = description,
        pedigree = pedigree,
        conductivity = conductivity,
    )
end

isdigitex(c::Char) = isdigit(c) || ((c >= '₀') && (c <= '₉'))

function remapdigits(str::AbstractString)
    res = str
    for repl in ('₀' + i => '0' + i for i = 0:9)
        res = replace(res, repl)
    end
    return res
end

"""
    Base.parse(
        ::Type{Material},
        expr::AbstractString;
        name::Union{AbstractString,Missing}=missing,
        properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
        atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
        density::Union{Missing, AbstractFloat}=missing,
        description::Union{Missing, AbstractString}=missing,
        pedigree::Union{Missing, AbstractString}=missing,
        conductivity::Union{Missing, Symbol}=missing, # :Conductor, :Semiconductor, :Insulator
    )::Material

Parse a Material from a string.
"""
function Base.parse(
    ::Type{Material},
    expr::AbstractString;
    name::Union{AbstractString,Missing} = missing,
    properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat} = missing,
    description::Union{Missing,AbstractString} = missing,
    pedigree::Union{Missing,AbstractString} = missing,
    conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
) where {V<:AbstractFloat}
    # Parses expressions like SiO2, Al2O3 or other simple (element qty) phrases
    function parseCompH1(expr::AbstractString)::Dict{PeriodicTable.Element,Int}
        parseSymbol(expr::AbstractString) =
            findfirst(z -> isequal(elements[z].symbol, expr), eachindex(elements))
        res = Dict{PeriodicTable.Element,Int}()
        start, idx = 1, collect(eachindex(expr))
        for i in eachindex(idx)
            if i < start
                continue
            elseif (i == start) || (i == start + 1) # Abbreviations are 1 or 2 letters
                if (i == start) && !isuppercase(expr[idx[i]]) # Abbrevs start with cap
                    error(
                        "Element abbreviations must start with a capital letter. $(expr[idx[i]])",
                    )
                end
                next = i + 1
                if (next > length(idx)) ||
                   isuppercase(expr[idx[next]]) ||
                   isdigitex(expr[idx[next]])
                    z = parseSymbol(expr[idx[start]:idx[i]])
                    if isnothing(z)
                        error("Unrecognized element parsing compound: $(expr[start:i])")
                    end
                    elm, cx = elements[z], 1
                    if (next <= length(idx)) && isdigitex(expr[idx[next]])
                        for stop = next:length(idx)
                            if (stop == length(idx)) || (!isdigitex(expr[idx[stop+1]]))
                                cx = parse(Int, remapdigits(expr[idx[next]:idx[stop]]))
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
    function parseCompH2(expr::AbstractString)::Dict{PeriodicTable.Element,Int}
        # println("Parsing: $(expr)")
        splt = split(expr, c -> (c in ('⋅', '.', '·')))
        if length(splt) > 1
            return mapreduce(parseCompH2, (a, b) -> merge(+, a, b), splt)
        end
        cx, start, stop = 0, -1, -1
        for i in eachindex(expr)
            if expr[i] == '('
                if cx == 0
                    start = i
                end
                cx += 1
            end
            if expr[i] == ')'
                cx -= 1
                if cx < 0
                    error("Unmatched right parenthesis.")
                elseif cx == 0
                    stop = i
                    break
                end
            end
        end
        if (start > 0) && (stop > start)
            res, q = parseCompH2(expr[nextind(expr, start):prevind(expr, stop)]), 1
            if (nextind(expr, stop) > lastindex(expr)) ||
               isdigitex(expr[nextind(expr, stop)])
                for i = stop:lastindex(expr)
                    if (nextind(expr, i) > lastindex(expr)) ||
                       (!isdigitex(expr[nextind(expr, i)]))
                        nxt = nextind(expr, stop)
                        q =
                            nxt <= lastindex(expr) ? parse(Int, remapdigits(expr[nxt:i])) :
                            1
                        stop = i
                        break
                    end
                end
            end
            for elm in keys(res)
                res[elm] *= q
            end
            if start > 1
                res = merge(+, res, parseCompH2(expr[1:prevind(expr, start)]))
            end
            if stop < lastindex(expr)
                res = merge(+, res, parseCompH2(expr[nextind(expr, stop):lastindex(expr)]))
            end
        else
            res = parseCompH1(expr)
        end
        # println("Parsing: $(expr) to $(res)")
        return res
    end
    function parseC(str::AbstractString)::Union{UncertainValue, Float64}
        str=strip(str)
        if startswith(str, "(") && endswith(str, ")")
            str=str[2:end-1]
        end
        uv = parse(UncertainValue, str)
        return σ(uv)==0.0 ? value(uv) : uv
    end

    # First split sums of Materials
    splt = split(expr, c -> c == '+')
    if length(splt) > 1
        return mapreduce(
            t -> parse(
                Material,
                strip(t),
                properties = properties,
                atomicweights = atomicweights,
                density = density,
                description = description,
                conductivity = conductivity,
                pedigree = pedigree,
            ),
            (a, b) -> sum(
                a,
                b,
                name = name,
                properties = properties,
                density = density,
                description = description,
                conductivity = conductivity,
                pedigree = pedigree,
            ),
            splt,
        )
    end
    # Second handle N*Material where N is a number
    p = findfirst(c -> (c == '*') || (c == '×') || (c == '⋅'), expr)
    if !isnothing(p)
        nm = ismissing(name) ? expr[nextind(expr, p):lastindex(expr)] : name
        return parseC(expr[firstindex(expr):prevind(expr, p)]) * parse(
            Material,
            expr[nextind(expr, p):lastindex(expr)],
            name = nm,
            properties = properties,
            atomicweights = atomicweights,
            density = density,
            description = description,
            conductivity = conductivity,
            pedigree = pedigree,
        )
    end
    # Then handle Material/N where N is a number
    p = findfirst(c -> c == '/', expr)
    if !isnothing(p)
        nm = ismissing(name) ? expr[firstindex(expr):previdx(expr, p)] : name
        return (1.0 / parse(Float64, expr[nextidx(expr, p):lastindex(expr)])) * parse(
            Material,
            expr[firstindex(expr):prevind(expr, p)],
            name = nm,
            properties = properties,
            atomicweights = atomicweights,
            density = density,
            description = description,
            conductivity = conductivity,
            pedigree = pedigree,
        )
    end
    # Finally parse material
    return atomicfraction(
        ismissing(name) ? expr : name,
        parseCompH2(expr);
        properties = properties,
        atomicweights = atomicweights,
        density = density,
        description = description,
        conductivity = conductivity,
        pedigree = pedigree,
    )
end

macro mat_str(str)
    parse(Material, str)
end

"""
    NeXLUncertainties.asa(::Type{DataFrame}, mat::Material)

Tabulate the composition of this Material as a DataFrame.  Columns for
material name, element abbreviation, atomic number, atomic weight, mass fraction,
normalized mass fraction, and atomic fraction. Rows for each element in mat.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, mat::Material)
    res = DataFrame(
        Material = Vector{String}(),
        Element = Vector{String}(),
        Z = Vector{Int}(),
        A = Vector{AbstractFloat}(),
        MassFrac = Vector{AbstractFloat}(),
        NormMass = Vector{AbstractFloat}(),
        AtomFrac = Vector{AbstractFloat}(),
    )
    af, tot = atomicfraction(mat), analyticaltotal(mat)
    for elm in sort(collect(keys(mat)))
        push!(
            res,
            (
                name(mat),
                symbol(elm),
                z(elm),
                a(elm, mat),
                mat[elm],
                mat[elm] / tot,
                af[elm],
            ),
        )
    end
    return res
end

"""
    NeXLUncertainties.asa(::Type{DataFrame}, mats::AbstractArray{Material}, mode=:MassFraction)

Tabulate the composition of a list of materials in a DataFrame.  One column
for each element in any of the materials.

    mode = :MassFraction | :NormalizedMassFraction | :AtomicFraction.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    mats::AbstractArray{<:Material},
    mode = :MassFraction,
)
    elms =
        length(mats) == 1 ? collect(keys(mats[1])) :
        Base.convert(Vector{Element}, sort(reduce(union, keys.(mats)))) # array of sorted Element
    cols = (Symbol("Material"), Symbol.(symbol.(elms))..., Symbol("Total")) # Column names
    empty = NamedTuple{cols}(
        map(c -> c == :Material ? Vector{String}() : Vector{AbstractFloat}(), cols),
    )
    res = DataFrame(empty) # Emtpy data frame with necessary columns
    for mat in mats
        vals = if mode == :AtomicFraction
            atomicfraction(mat)
        elseif mode == :NormalizedMassFraction
            normalizedmassfraction(mat)
        else
            massfraction(mat)
        end
        tmp = [
            name(mat),
            (value(get(vals, elm, 0.0)) for elm in elms)...,
            analyticaltotal(mat),
        ]
        push!(res, tmp)
    end
    return res
end

"""
    NeXLUncertainties.asa(::Type{LaTeXString}, mat::Material; parsename=true, order = :massfraction | :z)

Converts a `Material` into a `LaTeXString`.  `parsename` controls whether the material name is assumed to
be a parsable chemical formula (according to \\ce{...}).
"""
function NeXLUncertainties.asa(
    ::Type{LaTeXString},
    mat::Material;
    parsename = true,
    order = :massfraction,
)
    elms =
        order == massfraction ? #
        sort(collect(keys(mat)), lt = (e1, e2) -> mat[e1] > mat[e2]) :
        sort(collect(keys(mat)))
    cstr = join(
        ["\\ce{$(symbol(elm))}:\\num{$(round(mat[elm], digits=4))}" for elm in elms],
        ", ",
    )
    nm = parsename ? "\\ce{$(name(mat))}" : "\\mathrm{$(name(mat))}"
    return latexstring("$nm~:~\\left( $cstr \\mathrm{~by~mass} \\right)")
end

"""
    compare(unk::Material, known::Material)::DataFrame

Compare two compositions in a DataFrame.
"""
function compare(unk::Material, known::Material)::DataFrame
    um, km, z, kmf, rmf, dmf = Vector{String}(),
    Vector{String}(),
    Vector{String}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}()
    fmf, kaf, raf, daf, faf = Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}()
    afk, afr = atomicfraction(known), atomicfraction(unk)
    for elm in union(keys(known), keys(unk))
        push!(um, name(unk))
        push!(km, name(known))
        push!(z, elm.symbol)
        push!(kmf, known[elm])
        push!(rmf, value(unk[elm]))
        push!(dmf, value(known[elm]) - value(unk[elm]))
        push!(fmf, (value(known[elm]) - value(unk[elm])) / value(known[elm]))
        push!(kaf, value(get(afk, elm, 0.0)))
        push!(raf, value(get(afr, elm, 0.0)))
        push!(daf, value(get(afk, elm, 0.0)) - value(get(afr, elm, 0.0)))
        push!(
            faf,
            100.0 * (value(get(afk, elm, 0.0)) - value(get(afr, elm, 0.0))) /
            value(get(afk, elm, 0.0)),
        )
    end
    return DataFrame(
        Unkown = um,
        Known = km,
        Elm = z,
        Cknown = kmf,
        Cresult = rmf,
        ΔC = dmf,
        ΔCoC = fmf,
        Aknown = kaf,
        Aresult = raf,
        ΔA = daf,
        ΔAoA = faf,
    )
end

compare(unks::AbstractVector{<:Material}, known::Material) =
    mapreduce(unk -> compare(unk, known), append!, unks)


"""
    mac(mat::Material, xray::Union{Float64,CharXRay}, alg::Type{<:NeXLAlgorithm}=FFASTDB)::Float64

Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, xray::Union{Float64,CharXRay}, alg::Type{<:NeXLAlgorithm} = FFASTDB) =
    mapreduce(elm -> mac(elm, xray, alg) * value(mat[elm]), +, keys(mat))

function parsedtsa2comp(value::AbstractString)::Material
    try
        sp = split(value, ",")
        name = sp[1]
        mf, density = Dict{Element,Float64}(), missing
        for item in sp[2:end]
            if item[1] == '(' && item[end] == ')'
                sp2 = split(item[2:end-1], ":")
                mf[parse(Element, sp2[1])] = 0.01 * parse(Float64, sp2[2])
            else
                density = parse(Float64, item)
            end
        end
        return material(name, mf, density = density)
    catch err
        @warn "Error parsing composition $(value) - $(err)"
    end
end

function todtsa2comp(mat::Material)::String
    res = replace(name(mat), ',' => '_')
    elms = sort(collect(keys(mat)))
    for elm in elms
        qty = mat[elm]
        res *= ",($(symbol(elm)):$(round(100.0*qty,digits=2)))"
    end
    if haskey(mat.properties, :Density)
        res *= ",$(mat.properties[:Density])"
    end
    return res
end

"""
    compositionlibrary()::Dict{String, Material}

Load the internal compositon library.
"""
function compositionlibrary()::Dict{String,Material}
    result = Dict{String,Material}()
    path = dirname(pathof(@__MODULE__))
    df = CSV.File(joinpath(path, "..", "data", "composition.csv")) |> DataFrame
    for row in eachrow(df)
        name, density = row[1], row[2]
        elmc = collect(zip(element.(1:94), row[3:96])) # (i->getindex(row,i)).(3:96)))
        data = Dict{Element,Float64}(filter(a -> (!ismissing(a[2])) && (a[2] > 0.0), elmc))
        properties = Dict{Symbol,Any}()
        if !ismissing(density)
            properties[:Density] = density
        end
        m = material(name, data; properties = properties)
        result[name] = m
    end
    return result
end

"""
    z(mat::Material)


Computes the mean atomic number for a material. (Naive algorithm.)
"""
z(mat::Material) = sum(c * z for (z, c) in mat.massfraction)


"""
    z(mat::Material)


Computes the mean atomic weight for a material. (Naive algorithm.)
"""
a(mat::Material) =
    sum(haskey(mat.a, z) ? mat.a[z] : c * a(elements[z]) for (z, c) in mat.massfraction)
    

"""
    Statistics.mean(mats::AbstractArray{<:Material})

If the mass fractions for all the elements in all `mats` have non-zero uncertainties
then the variance weighted mean is calculated and the result will have associated 
uncertainties.  Otherwise, the straight floating-point mean is calculated and
the result won't have uncertainties.  This is because even a single value with
zero uncertainty will poison the variance weighted mean (produce a NaN).
"""
function Statistics.mean(mats::AbstractArray{<:Material})
    zs = mapreduce(m->elms(m), union!, mats)
    nm = if length(mats)>5
        "mean[$(name(mats[1])) + $(length(mats)-1) others]"
    else
        "mean[$(join(name.(mats), ", "))]"
    end
    function mm(z)
        return if all(σ(mat[z]) > 0.0 for mat in mats)
            mean(UncertainValue[uv(mat[z]) for mat in mats])
        else
            mean(Float64[value(mat[z]) for mat in mats])
        end
    end
    return material(nm, (z=>mm(z) for z in zs)...)
end

"""
    Base.rand(::Type{Material}, zs::AbstractUnitRange{Int}=1:20)::Material

Generate a randomize material.
"""
function Base.rand(::Type{Material}, zs::AbstractUnitRange{Int}=1:20)::Material{UncertainValue,Float64}
    sum, mfs = 0.0, Dict{Element,UncertainValue}()
    while sum<1.0
        z, r=Base.rand(zs), Base.rand()
        if !haskey(mfs, elements[z])
            v = min(r, 1.0-sum)
            mfs[elements[z]] = uv(v,Base.rand()*0.1*v)
            sum+=r
        end
    end
    return material("random", mfs)
end

"""
    Base.similar(mat::Material{UncertainValue, <:AbstractFloat}, n::Integer)::Vector{Material{UncertainValue,Float64}}

Generate `n` Materials similar to `mat` using the uncertainties in `mat` as 
your guide of dispersion.  The mass-fractions of `mat` must
be defined as `UncertainValue`s.
"""
function Base.similar(mat::Material{UncertainValue, <:AbstractFloat}, n::Integer)::Vector{Material{UncertainValue,Float64}}
    return [ material(
        "Like[$(name(mat)), $i]", 
        Dict(z=>uv(value(mat[z])+(1.0-2.0*Base.rand())*σ(mat[z]), (0.9+0.2*Base.rand()*σ(mat[z]))) 
            for z in keys(mat))
    ) for i in 1:n ]
end


