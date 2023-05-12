# Defines the Material struct and functions on it.
using Unitful
using Printf
using DataFrames
using LaTeXStrings
import Base.rand
import Statistics

"""
Holds basic data about a material including name, composition in mass fraction and optional propreties.

By default, Material assumes nominal terrestrial atomic weights.  However, it is possible to assign custom
atomic weights on a per element-basis for non-terrestrial materials.

The mass fraction and atomic weight are immutable but the `Properties` can be modified.

    Material(
        name::AbstractString,
        massfrac::AbstractDict{Element,U},
        atomicweights::AbstractDict{Element,V} = Dict{Element,Float64}(),
        properties::AbstractDict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where { U <: AbstractFloat, V <: AbstractFloat }


**Properties**

    :Density # Density in g/cm³
    :Description # Human friendly
    :Pedigree #  Quality indicator for compositional data ("SRM-XXX", "CRM-XXX", "NIST K-Glass", "Stoichiometry", "Wet-chemistry by ???", "WDS by ???", "???")
    :Conductivity = :Insulator | :Semiconductor | :Conductor
    :OtherUserProperties # Other properties can be defined as needed
"""
struct Material{U<:AbstractFloat,V<:AbstractFloat}
    name::String
    massfraction::Dict{Element,U}
    a::Dict{Element,V} # Optional: custom atomic weights for the keys in this Material
    properties::Dict{Symbol,Any} # :Density, :Description, :Pedigree, :Conductivity, ... + user defined

    """
        Material(
            name::AbstractString,
            massfrac::AbstractDict{Element,U},
            atomicweights::AbstractDict{Element,V} = Dict{Element,Float64}(),
            properties::AbstractDict{Symbol,Any} = Dict{Symbol,Any}(),
        ) where { U <: AbstractFloat, V <: AbstractFloat }
    """
    function Material(
        name::AbstractString,
        massfrac::AbstractDict{Element,U},
        atomicweights::AbstractDict{Element,V}=Dict{Element,Float64}(),
        properties::AbstractDict{Symbol,Any}=Dict{Symbol,Any}(),
    ) where {U<:AbstractFloat,V<:AbstractFloat}
        if sum(value.(values(massfrac))) > 50.0
            @warn "The sum mass fraction is $(sum(values(massfrac))) which is much larger than unity."
        end
        new{U,V}(name, massfrac, atomicweights, properties)
    end
end

const NULL_MATERIAL = Material("Null Material", Dict{Element,Float64}())
"""
    rename(mat::Material, newname::AbstractString)

Creates a replica of `mat` but with a new name.
"""
rename(mat::Material, newname::AbstractString) = Material(newname, mat.massfraction, mat.a, mat.properties)

Base.copy(m::Material) =
    Material(m.name, copy(m.massfraction), copy(m.a), copy(m.properties))

properties(mat::Material) = mat.properties

"""
    elms(mat::Material)

The elements with mass fraction ≠ 0.0 in `mat`.
"""
elms(mat::Material) = keys(mat.massfraction)

"""
    ispure(mat::Material)

Does `mat` represent a single element.
"""
ispure(mat::Material) = length(mat.massfraction) == 1

function Base.:*(k::AbstractFloat, mat::Material)::Material
    mf = Dict(el => q * k for (el, q) in mat.massfraction)
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
    name::Union{AbstractString,Missing}=missing,
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing # :Conductor, :Semiconductor, :Insulator
)::Material
    plus(v1::AbstractFloat, v2::AbstractFloat) =
        (σ(v1) == 0.0) && (σ(v2) == 0.0) ? value(v1) + value(v2) : uv(value(v1) + value(v2), sqrt(σ(v1)^2 + σ(v2)^2))
    mf = Dict{Element,AbstractFloat}(elm => plus(mat1[elm], mat2[elm]) for elm in union(elms(mat1), elms(mat2)))
    aw = Dict{Element,Float64}()
    for elm in union(keys(mat1.a), keys(mat2.a))
        aw[elm] =
            (value(mat1[elm]) + value(mat2[elm])) / (value(mat1[elm]) / a(elm, mat1) + value(mat2[elm]) / a(elm, mat2))
    end
    name = ismissing(name) ? "$(mat1.name)+$(mat2.name)" : name
    return material(
        name,
        mf;
        properties=properties,
        atomicweights=aw,
        density=density,
        description=description,
        conductivity=conductivity,
        pedigree=pedigree
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
    data::Dict{Material,<:AbstractFloat};
    name::Union{AbstractString,Missing}=missing,
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing # :Conductor, :Semiconductor, :Insulator
)::Material
    assign(val, prop, props) =
        if !ismissing(val)
            props[prop] = val
        end
    res = mapreduce((a, b) -> sum(a, b), data) do (mat, f)
        f * mat
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
    atoms_per_g(mat::Material)
    atoms_per_g(mat::Material, elm::Element)

Compute the number of atoms of `elm` in 1 gram of `mat`.
"""
atoms_per_g(elm::Element) = ustrip(NoUnits, AvogadroConstant / (a(elm) * u"1/mol"))
atoms_per_g(mat::Material, elm::Element) = ustrip(NoUnits, mat[elm] * AvogadroConstant / (a(elm, mat) * u"1/mol"))
atoms_per_g(mat::Material) = sum(elm -> atoms_per_g(mat, elm), keys(mat))

"""
    material(
        name::AbstractString,
        massfrac::Dict{Element,U};
        properties::Union{Missing,Dict{Symbol,Any}} = missing,
        atomicweights::Union{Missing, Dict{Element,Float64}} = missing,
        density::Union{Missing,AbstractFloat} = missing,
        description::Union{Missing,AbstractString} = missing,
        pedigree::Union{Missing,AbstractString} = missing,
        conductivity::Union{Missing,Symbol} = missing, # :Conductor, :Semiconductor, :Insulator
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
    properties::Union{Missing,Dict{Symbol,Any}}=missing,
    atomicweights::Union{Missing,Dict{Element,Float64}}=missing,
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing # :Conductor, :Semiconductor, :Insulator
) where {U<:AbstractFloat}
    props = ismissing(properties) ? Dict{Symbol,Any}() : copy(properties)
    atomicweights = ismissing(atomicweights) ? Dict{Element,Float64}() : copy(atomicweights)
    (!ismissing(density)) && ((props[:Density] = density) == density)
    (!ismissing(description)) && ((props[:Description] = description) == description)
    (!ismissing(pedigree)) && ((props[:Pedigree] = pedigree) == pedigree)
    (!ismissing(conductivity)) && ((props[:Conductivity] = conductivity) == conductivity)
    return Material(name, massfrac, atomicweights, props)
end

material(
    name::AbstractString,
    massfrac::Pair{Element,U}...;
    properties::Union{Missing,Dict{Symbol,Any}}=missing,
    atomicweights::Union{Missing,Dict{Element,Float64}}=missing,
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing # :Conductor, :Semiconductor, :Insulator
) where {U<:AbstractFloat} = # 
    material(
        name,
        Dict(massfrac);
        properties=properties,
        atomicweights=atomicweights,
        density=density,
        description=description,
        pedigree=pedigree,
        conductivity=conductivity
    )

"""
     material(str::String, density::Float64)
     nargs...
    pure(elm::Element)

Construct a Material to represent a pure element.

Example:

    > pure(n"Fe")
"""
pure(elm::Element) =
    material("Pure $(symbol(elm))", Dict{}(elm => 1.0), density=density(elm))

function Base.show(io::IO, mat::Material)
    res = "$(name(mat))["
    res *= join(
        (
            @sprintf("%s=%0.4f", el.symbol, value(mat[el])) for el in sort(collect(keys(mat.massfraction)))
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
    return material(name, tmp, density=density)
end

Base.convert(::Type{Material{Float64,Float64}}, comp::Material{Float64,Float64}) = comp
function Base.convert(::Type{Material{Float64,Float64}}, comp::Material)::Material{Float64,Float64}
    mfs = Dict(el => Float64(value(mf)) for (el, mf) in comp.massfraction)
    aws = Dict(el => Float64(value(aa)) for (el, aa) in comp.a)
    return Material(name(comp), mfs, aws, comp.properties)
end
Base.convert(::Type{Material{Float64,Float64}}, ::Nothing) = nothing


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

"""
    nonneg(mat::Material{U,V}, elm::Element)::U where {U<:AbstractFloat,V<:AbstractFloat}
    nonneg(mat::Material{UncertainValue,V}, elm::Element)::Float64 where {V<:AbstractFloat}
    nonneg(mat::Material)::Material

Returns the mass fraction of `elm::Element` truncated to be non-negative.  Negative values
are returned as 0.0. Positive values are returned as is.
"""
function nonneg(mat::Material{U,V}, elm::Element)::U where {U<:AbstractFloat,V<:AbstractFloat}
    max(zero(U), value(mat[elm]))
end
function nonneg(mat::Material{UncertainValue,V}, elm::Element)::Float64 where {V<:AbstractFloat}
    max(0.0, value(mat[elm]))
end

nonneg(mat::Material) = #
    Material(mat.name, Dict(el => nonneg(mat, el) for el in keys(mat.massfraction)), mat.a, mat.properties)

"""
    normalizedmassfraction(mat::Material)::Dict{Element, AbstractFloat}

Return the normalized mass fraction as a Dict{Element, AbstractFloat}.  Negative values
are set to zero.
"""
function normalizedmassfraction(mat::Material)::Dict{Element,AbstractFloat}
    n = analyticaltotal(mat)
    return Dict(elm => nonneg(mat, elm) / n for elm in keys(mat))
end


"""
    normalized(mat::Material{U,V}, elm::Element)

Returns the mass fraction of 'elm::Element' such that the returned value is non-negative
and the sum of all values is unity.
"""
function normalized(mat::Material{U,V}, elm::Element) where {U<:AbstractFloat,V<:AbstractFloat}
    nonneg(mat, elm) / analyticaltotal(mat)
end

"""
    asnormalized(mat::Material, n=1.0)::Material

Convert the Material to a normalized Material form.  Negative mass fractions
are set to zero before normalization.
"""
function asnormalized(mat::Material{U,V}, n=one(U)) where {U<:AbstractFloat,V<:AbstractFloat}
    at = analyticaltotal(mat)
    if isapprox(at, n, rtol=1.0e-8) && startswith(mat.name, "N[")
        return mat
    else
        if at == 0.0
            at = 1.0
        end
        return Material(
            "N[$(name(mat)),$(n)]",
            Dict{Element,U}(elm => n * (nonneg(mat, elm) / at) for elm in keys(mat)),
            mat.a,
            copy(mat.properties),
        )
    end
end

"""
    Base.isapprox(mat1::Material, mat2::Material; atol = 1.0e-4)

Are these Material(s) equivalent to within `atol`?
"""
function Base.isapprox(mat1::Material, mat2::Material; atol=1.0e-4)
    return all(
        isapprox(value(mat1[elm]), value(mat2[elm]), atol=atol) && #
        isapprox(σ(mat1[elm]), σ(mat2[elm]), atol=atol)
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
    atomicfraction(mat::Material{U,V})::Dict{Element,U}

Return the composition in atomic fraction representation.
"""
function atomicfraction(mat::Material{U,V})::Dict{Element,U} where {U<:AbstractFloat,V<:AbstractFloat}
    norm = sum(value(mf) / value(a(elm, mat)) for (elm, mf) in mat.massfraction; init=zero(Float64))
    return Dict(elm => (mf / value(a(elm, mat))) / norm for (elm, mf) in mat.massfraction)
end

"""
    analyticaltotal(mat::Material)

Return the sum of the positive mass fractions.
"""
function analyticaltotal(mat::Material{U,V})::U where {U<:AbstractFloat,V<:AbstractFloat}
    sum(val -> value(val) < 0.0 ? zero(U) : val, values(mat.massfraction); init=zero(U))
end

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
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V}=Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing # :Conductor, :Semiconductor, :Insulator
) where {U<:Real,V<:AbstractFloat} = atomicfraction(
    name,
    Dict(atomfracs);
    properties=properties,
    atomicweights=atomicweights,
    density=density,
    description=description,
    pedigree=pedigree,
    conductivity=conductivity
)

function atomicfraction(
    name::AbstractString,
    atomfracs::Dict{Element,U};
    properties::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    atomicweights::Dict{Element,V}=Dict{Element,Float64}(),
    density::Union{Missing,AbstractFloat}=missing,
    description::Union{Missing,AbstractString}=missing,
    pedigree::Union{Missing,AbstractString}=missing,
    conductivity::Union{Missing,Symbol}=missing # :Conductor, :Semiconductor, :Insulator
) where {U<:Real,V<:AbstractFloat}
    aw(elm) = get(atomicweights, elm, a(elm))
    norm = sum(af * aw(elm) for (elm, af) in atomfracs)
    massfracs = Dict(elm => (af * aw(elm) / norm) for (elm, af) in atomfracs)
    return material(
        name,
        massfracs;
        atomicweights=atomicweights,
        properties=properties,
        density=density,
        description=description,
        pedigree=pedigree,
        conductivity=conductivity
    )
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
        "Material" => [name(mat) for _ in els],
        "Element" => [el.symbol for el in els],
        "Z" => [z(el) for el in els],
        "A" => [a(el, mat) for el in els],
        "C(z)" => [mat[el] for el in els],
        "Norm[C(z)]" => [nmf[el] for el in els],
        "A(z)" => [af[el] for el in els]
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
    mode=:MassFraction,
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
    parsename=true,
    order=:massfraction
)
    elms = if order == :massfraction
        sort(collect(keys(mat)), lt=(e1, e2) -> mat[e1] > mat[e2])
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
    els = collect(union(keys(known), keys(unk)))
    return DataFrame(
        Symbol("Material 1") => [name(unk) for _ in els],
        Symbol("Material 2") => [name(known) for _ in els],
        Symbol("Elm") => [symbol(el) for el in els],
        Symbol("C₁(z)") => [known[el] for el in els],
        Symbol("C₂(z)") => [value(unk[el]) for el in els],
        Symbol("ΔC") => [value(known[el]) - value(unk[el]) for el in els],
        Symbol("ΔC/C") => map(els) do el
            (value(known[el]) - value(unk[el])) / #
            max(value(known[el]), value(unk[el]))
        end,
        Symbol("A₁(z)") => [value(get(afk, el, 0.0)) for el in els],
        Symbol("A₂(z)") => [value(get(afr, el, 0.0)) for el in els],
        Symbol("ΔA") => [value(get(afk, el, 0.0)) - value(get(afr, el, 0.0)) for el in els],
        Symbol("ΔA/A") => map(els) do el
            (value(get(afk, el, 0.0)) - value(get(afr, el, 0.0))) / #
            max(value(get(afk, el, 0.0)), value(get(afr, el, 0.0)))
        end, copycols=false
    )
end

compare(unks::AbstractVector{<:Material}, known::Material) =
    mapreduce(unk -> compare(unk, known), append!, unks)


"""
    mac(mat::Material, xray::Union{Float64,CharXRay}, alg::Type{<:NeXLAlgorithm}=DefaultAlgorithm)::Float64

Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, energy::Float64, alg::Type{<:NeXLAlgorithm}=DefaultAlgorithm) =
    sum(mat.massfraction) do (elm, mf)
        mac(elm, energy, alg) * max(0.0, value(mf))
    end
mac(mat::Material, xray::CharXRay, alg::Type{<:NeXLAlgorithm}=DefaultAlgorithm) =
    sum(mat.massfraction) do (elm, mf)
        mac(elm, xray, alg) * max(0.0, value(mf))
    end

listcustommacs(mat::Material) = listcustommacs(keys(mat))

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
        return material(name, mf, density=density)
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
    z(mat::Material) = z(Donovan2002, mat)
    z(::Type{NaiveZ|ElectronFraction|AtomicFraction}, mat::Material)
    z(::Type{ElasticFraction}, mat::Material, e::AbstractFloat)
    z(::Type{Donovan2002}, mat::Material; exponent=0.667)

Compute the mean atomic number for a material.

    Algorithms:
      * NaiveZ - Mass fraction averaging
      * AtomFraction - Atom fraction averaging
      * ElectronFraction - Simple electron fraction averaging
      * ElasticFraction - Scattering cross-section averaged
      * Donovan2002 - Yukawa/Donovan modified exponent electron fraction averaging

For more details see Mean Z algorithm in J.J. Donovan, N.E. Pingitore, Microsc. Microanal. 2002 ; 8 , 429
(also see Microsc. Microanal. 27 (Suppl 1), 2021))
"""
z(mat::Material) = z(Donovan2002, mat)

struct NaiveZ <: NeXLAlgorithm end
z(::Type{NaiveZ}, mat::Material) = sum(c * z(elm) for (elm, c) in mat.massfraction)

"""
Donovan's recommended material Z model

For more details see Mean Z algorithm in J.J. Donovan, N.E. Pingitore, Microsc. Microanal. 2002 ; 8 , 429
(also see Microsc. Microanal. 27 (Suppl 1), 2021))
"""
struct Donovan2002 <: NeXLAlgorithm end
function z(::Type{Donovan2002}, mat::Material; exponent=0.667)
    af = atomicfraction(mat)
    return sum(NeXLUncertainties.value(a) * z(elm)^(1.0 + exponent) for (elm, a) in af) / #
           sum(NeXLUncertainties.value(a) * z(elm)^exponent for (elm, a) in af)
end


"""
Electronic fraction material Z model

For more details see Mean Z algorithm in J.J. Donovan, N.E. Pingitore, Microsc. Microanal. 2002 ; 8 , 429
(also see Microsc. Microanal. 27 (Suppl 1), 2021))
"""
struct ElectronFraction <: NeXLAlgorithm end

function z(::Type{ElectronFraction}, mat::Material)
    af = atomicfraction(mat)
    ef(elm) = af[elm] * z(elm) / sum(el2 -> af[el2] * z(el2), keys(mat)) # Donovan2002 Eq 3
    return sum(elm -> z(elm) * ef(elm), keys(mat))
end

"""
Elastic fraction material Z model

For more details see Mean Z algorithm in J.J. Donovan, N.E. Pingitore, Microsc. Microanal. 2002 ; 8 , 429
(also see Microsc. Microanal. 27 (Suppl 1), 2021))
"""
struct ElasticFraction <: NeXLAlgorithm end

function z(::Type{ElasticFraction}, mat::Material, e::AbstractFloat)
    function σE(Z, E) # E in keV
        α = 3.4e-3 * Z^0.67 / E
        return 5.21e-21 * (Z / E)^2 * (4π) / (α * (1.0 + α)) * ((E + 0.511e3) / (E + 2.0 * 0.511e3))^2
    end
    af = atomicfraction(mat)
    σf(elm) = af[elm] * σE(z(elm), 0.001 * e) / sum(el -> af[el] * σE(z(el), 0.001 * e), keys(mat))
    return sum(elm -> z(elm) * σf(elm), keys(mat))
end

"""
Naive atomic fraction material Z model

For more details see Mean Z algorithm in J.J. Donovan, N.E. Pingitore, Microsc. Microanal. 2002 ; 8 , 429
(also see Microsc. Microanal. 27 (Suppl 1), 2021))
"""
struct AtomicFraction <: NeXLAlgorithm end

function z(::Type{AtomicFraction}, mat::Material)
    af = atomicfraction(mat)
    sum(elm -> af[elm] * z(elm), keys(mat))
end


"""
    a(mat::Material)

Computes the mean atomic weight for a material.
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
    els = mapreduce(m -> keys(m), union!, mats, init=Set{Element}())
    nm = if length(mats) > 5
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
    return material(nm, Dict(el => mm(el) for el in els))
end

"""
    Base.rand(::Type{Material}, zs::AbstractUnitRange{Int}=1:20)::Material

Generate a randomize material.
"""
function Base.rand(::Type{Material}, zs::AbstractUnitRange{Int}=1:20)::Material{UncertainValue,Float64}
    sum, mfs = 0.0, Dict{Element,UncertainValue}()
    while sum < 1.0
        z, r = Base.rand(zs), Base.rand()
        if !haskey(mfs, elements[z])
            v = min(r, 1.0 - sum)
            mfs[elements[z]] = uv(v, Base.rand() * 0.1 * v)
            sum += r
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
function Base.similar(mat::Material{UncertainValue,<:AbstractFloat}, n::Integer)::Vector{Material{UncertainValue,Float64}}
    return [material(
        "Like[$(name(mat)), $i]",
        Dict(el => uv(value(mat[el]) + (1.0 - 2.0 * Base.rand()) * σ(mat[el]), (0.9 + 0.2 * Base.rand() * σ(mat[el])))
             for el in keys(mat))
    ) for i in 1:n]
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


"""
    NeXLUncertainties.asa(Dict, mat::Material)

Convert a `Material` into a `Dict{String, Any}` as is suitable
for conversion to JSON.
"""
function NeXLUncertainties.asa(::Type{Dict}, mat::Material)
    res = Dict{String,Any}()
    res["Name"] = name(mat)
    res["MassFraction"] = Dict(symbol(elm) => value(mat[elm]) for elm in keys(mat))
    if !isempty(mat.a)
        res["AtomicWeight"] = Dict(symbol(elm) => mat[elm] for elm in keys(mat.a))
    end
    if haskey(mat.properties, :Density)
        res["Density"] = mat.properties[:Density]
    end
    if haskey(mat.properties, :Description)
        res["Description"] = mat.properties[:Description]
    end
    # AtomicFraction and NormalizedMassFraction are for seach purposes only...
    nmf = normalizedmassfraction(mat)
    res["NormalizedMassFraction"] = Dict(symbol(elm) => nmf[elm] for elm in keys(mat))
    af = atomicfraction(mat)
    res["AtomFraction"] = Dict(symbol(elm) => value(get(af, elm, zero(valtype(af)))) for elm in keys(mat))
    return res
end

"""
    Material(d::Dict)

Construct a `Material` from a `Dict` created by `NeXLUncertainties.asa(Dict, mat::Material)`
"""
function Material(d::Dict{String,Any})
    massfrac = Dict(
        parse(Element, elm) => q for (elm, q) in d["MassFraction"]
    )
    aw = haskey(d, "AtomicWeight") ? Dict(
        parse(Element, elm) => q for (elm, q) in d["AtomicWeight"]
    ) : Dict{Element,Float64}()
    props = Dict{Symbol,Any}()
    if haskey(d, "Density")
        props[:Density] = d["Density"]
    end
    if haskey(d, "Description")
        props[:Description] = d["Description"]
    end
    Material(d["Name"], massfrac, aw, props)
end
