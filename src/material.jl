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
    massfraction::Dict{Element,U}
    a::Dict{Element,V} # Optional: custom atomic weights for the keys in this Material
    
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
        massfrac::Dict{Element,U},
        atomicweights::Dict{Element,V} = Dict{Element,Float64}(),
        properties::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where {U<:AbstractFloat, V<:AbstractFloat}
        if sum(value.(values(massfrac))) > 10.0
            @warn "The sum mass fraction is $(sum(values(massfrac))) which is much larger than unity."
        end
        new{U,V}(name, properties, massfrac, atomicweights)
    end
end

const NULL_MATERIAL = Material("Null Material",Dict{Element,Float64}())
"""
    rename(mat::Material, newname::AbstractString)

Creates a replica of `mat` but with a new name.
"""
rename(mat::Material, newname::AbstractString) = Material(newname, mat.massfraction, mat.a, mat.properties)

Base.copy(m::Material) =
    Material(m.name, copy(m.massfraction), copy(m.a), copy(m.properties))

"""
    elms(mat::Material)

The elements with mass fraction ≠ 0.0 in `mat`.
"""
elms(mat::Material) = keys(mat.massfraction)

"""
    ispure(mat::Material)

Does `mat` represent a single element.
"""
ispure(mat::Material) = length(mat.massfraction)==1

function Base.:*(k::AbstractFloat, mat::Material)::Material
    mf = Dict( el => q * k for (el, q) in mat.massfraction )
    return Material("$(k)⋅$(mat.name)", mf, copy(mat.a), copy(mat.properties))
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
mixtures of compounds.  Ultimately, expressions like `mat"0.5*Al2O3+0.5*MgO"` or equivalently
`0.5*mat"Al2O3"+0.5*mat"MgO"` are computed using `sum(...)`.
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
    mf = Dict{Element,AbstractFloat}( elm => plus(mat1[elm], mat2[elm]) for elm in union(elms(mat1), elms(mat2)))
    aw = Dict{Element,Float64}()
    for elm in union(keys(mat1.a), keys(mat2.a))
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
density(mat::Material) = mat[:Density]

"""
    description(mat::Materail)

The :Description property.
"""
description(mat::Material) = mat[:Description]

"""
    pedigree(mat::Material)

The :Pedigree property.
"""
pedigree(mat::Material) = mat[Pedigree]

"""
    atoms_per_cm³(mat::Material, elm::Element) =

Number of atoms per cm³ of the specified Element in the specified Material.  The Material must define
the `:Density` property.

    atoms_per_cm³(mat::Material)

Total number of atoms per cm³ for all elements in mat. 
"""
atoms_per_cm³(mat::Material, elm::Element) = mat[:Density] * atoms_per_g(mat, elm)
atoms_per_cm³(mat::Material) = sum(atoms_per_cm³(mat, elm) for elm in keys(mat))


"""
    atoms_per_g(elm::Element)
    atoms_per_g(mat::Material, elm::Element)

Compute the number of atoms of `elm` in 1 gram of `mat`.
"""
atoms_per_g(elm::Element) =  ustrip(NoUnits, AvogadroConstant / (a(elm)*u"1/mol"))
atoms_per_g(mat::Material, elm::Element) = ustrip(NoUnits, mat[elm] * AvogadroConstant / (a(elm, mat)*u"1/mol"))

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
    return Material(name, massfrac, atomicweights, props)
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
            @sprintf("%s=%0.4f", elm.symbol, value(mf)) for (elm, mf) in mat.massfraction
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
a(elm::Element, mat::Material) = get(mat.a, elm, a(elm))

Base.getindex(mat::Material{U,V}, elm::Element) where {U<:AbstractFloat,V<:AbstractFloat} =
    get(mat.massfraction, elm, zero(U))
Base.getindex(mat::Material{U,V}, z::Int) where {U<:AbstractFloat,V<:AbstractFloat} =
    get(mat.massfraction, elements[z], zero(U))

Base.getindex(mat::Material, sym::Symbol) = getindex(mat.properties, sym)
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
    return Dict( elm => nonneg(mat, elm) / n for elm in keys(mat))
end

normalized(mat::Material, elm::Element) = nonneg(mat, elm) / analyticaltotal(mat)

"""
    asnormalized(mat::Material, n=1.0)::Material

Convert the Material to a normalized Material form.  Negative mass fractions
are set to zero before normalization.
"""
function asnormalized(mat::Material, n = 1.0)
    if !isempty(mat.massfraction)
        at = analyticaltotal(mat)
        if isapprox(at, n, rtol = 1.0e-8) && startswith(name(mat), "N[")
            return mat
        else
            return Material(
                "N[$(name(mat)),$(n)]",
                Dict( elm => n * nonneg(mat, elm) / max(1.0e-8, at) for elm in keys(mat)),
                mat.a,
                copy(mat.properties),
            )
        end
    else
        return mat
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
massfraction(mat::Material)::Dict{Element,AbstractFloat} = copy(mat.massfraction)

"""
    Base.keys(mat::Material)

Return an interator over the elements in the Material.
"""
Base.keys(mat::Material) = keys(mat.massfraction)

"""
    labeled(mat::Material)

Transform the mass fraction representation of a material into a Dict{MassFractionLabel,AbstractFloat}"""
labeled(mat::Material) =
    Dict(MassFractionLabel(name(mat), elm) => mf for (elm, mf) in mat.massfraction)

"""
    atomicfraction(mat::Material)::Dict{Element,AbstractFloat}

Return the composition in atomic fraction representation.
"""
function atomicfraction(mat::Material)::Dict{Element,AbstractFloat}
    norm = sum(mf / a(elm, mat) for (elm, mf) in mat.massfraction)
    return Dict( elm => (mf / a(elm, mat)) / norm for (elm, mf) in mat.massfraction )
end

"""
    analyticaltotal(mat::Material)

Return the sum of the positive mass fractions.
"""
analyticaltotal(mat::Material) = sum(elm->nonneg(mat, elm), keys(mat))

"""
    haskey(mat::Material, elm::Element)
    haskey(mat::Material, z::Int)

Does this material contain this element?
"""
Base.haskey(mat::Material, elm::Element) = haskey(mat.massfraction, elm)
Base.haskey(mat::Material, z::Int) = haskey(mat.massfraction, elements[z])

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
    massfracs = Dict(elm => (aw(elm) / norm) * af for (elm, af) in atomfracs)
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
        res = Dict{Element,Int}()
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
                    elm = get(elements, Symbol(expr[idx[start]:idx[i]]), nothing)
                    if isnothing(elm)
                        error("Unrecognized element parsing compound: $(expr[start:i])")
                    end
                    cx = 1
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
    af, nmf = atomicfraction(mat), normalizedmassfraction(mat)
    els = sort(collect(keys(mat)))
    return DataFrame(
        "Material" => [ name(mat) for _ in els ],
        "Element" => [ el.symbol for el in els ],
        "Z" => [ z(el) for el in els ],
        "A" => [ a(el, mat) for el in els ],
        "C(z)" => [ mat[el] for el in els ],
        "Norm[C(z)]" => [ nmf[el] for el in els ],
        "A(z)" => [ af[el] for el in els ]
    )
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
    elms = sort(collect(union(keys.(mats)...))) # array of sorted Element
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
    elms = if order == :massfraction
        sort(collect(keys(mat)), lt = (e1, e2) -> mat[e1] > mat[e2])
    else
        sort(collect(keys(mat)))
    end
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
    afk, afr = atomicfraction(known), atomicfraction(unk)
    els = collect(union(keys(known), keys(unk) ))
    return DataFrame(
        Symbol("Material 1") => [ name(unk) for _ in els ],
        Symbol("Material 2") => [ name(known) for _ in els ],
        Symbol("Elm") => [ symbol(el) for el in els ],
        Symbol("C₁(z)") => [ known[el] for el in els ],
        Symbol("C₂(z)") => [ value(unk[el]) for el in els ],
        Symbol("ΔC") => [ value(known[el]) - value(unk[el]) for el in els ],
        Symbol("ΔC/C") => map(els) do el
             (value(known[el]) - value(unk[el])) / #
              max(value(known[el]),value(unk[el]))
        end,
        Symbol("A₁(z)") => [ value(get(afk, el, 0.0)) for el in els ],
        Symbol("A₂(z)") => [ value(get(afr, el, 0.0)) for el in els ],
        Symbol("ΔA") => [ value(get(afk, el, 0.0)) - value(get(afr, el, 0.0)) for el in els ],
        Symbol("ΔA/A") => map(els)  do el
            (value(get(afk, el, 0.0)) - value(get(afr, el, 0.0))) / #
             max(value(get(afk, el, 0.0)),value(get(afr, el, 0.0))) 
        end
    )
end

compare(unks::AbstractVector{<:Material}, known::Material) =
    mapreduce(unk -> compare(unk, known), append!, unks)


"""
    mac(mat::Material, xray::Union{Float64,CharXRay}, alg::Type{<:NeXLAlgorithm}=FFASTDB)::Float64

Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, energy::Float64, alg::Type{<:NeXLAlgorithm} = FFASTDB) =
    sum(zc->mac(zc[1], energy, alg) * value(zc[2]), mat.massfraction) 
mac(mat::Material, xray::CharXRay, alg::Type{<:NeXLAlgorithm} = FFASTDB) =
    mac(mat, energy(xray), alg)

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
z(mat::Material) = sum(c * z(elm) for (elm, c) in mat.massfraction)


"""
    a(mat::Material)

Computes the mean atomic weight for a material. (Naive algorithm.)
"""
a(mat::Material) = sum(c * a(elm, mat) for (elm, c) in mat.massfraction)
    

"""
    Statistics.mean(mats::AbstractArray{<:Material})

If the mass fractions for all the elements in all `mats` have non-zero uncertainties
then the variance weighted mean is calculated and the result will have associated 
uncertainties.  Otherwise, the straight floating-point mean is calculated and
the result won't have uncertainties.  This is because even a single value with
zero uncertainty will poison the variance weighted mean (produce a NaN).
"""
function Statistics.mean(mats::AbstractArray{<:Material})
    els = mapreduce(m->keys(m), union!, mats, init=Set{Element}())
    nm = if length(mats)>5
        "mean[$(name(mats[1])) + $(length(mats)-1) others]"
    else
        "mean[$(join(name.(mats), ", "))]"
    end
    function mm(el)
        return if all(σ(mat[el]) > 0.0 for mat in mats)
            mean(UncertainValue[uv(mat[el]) for mat in mats])
        else
            mean(Float64[value(mat[el]) for mat in mats])
        end
    end
    return material(nm, Dict(el=>mm(el) for el in els))
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
        Dict(el=>uv(value(mat[el])+(1.0-2.0*Base.rand())*σ(mat[el]), (0.9+0.2*Base.rand()*σ(mat[el]))) 
            for el in keys(mat))
    ) for i in 1:n ]
end

"""
    delete(mat::Material, elm::Element)::Material
    delete(mat::Material, elm::AbstractVector{Element})::Material

Constructs a new Material from `mat` with `elm` removed.
"""
function delete(mat::Material, elm::Element)
    res = copy(mat)
    delete!(res.massfraction, elm)
    delete!(res.a, elm)
    return res
end
function delete(mat::Material, els::AbstractVector{Element})
    res = copy(mat)
    for elm in els
        delete!(res.massfraction, elm)
        delete!(res.a, elm)
    end
    return res
end
