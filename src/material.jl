# Defines the Material struct and functions on it.
using Unitful
using PeriodicTable
using Printf

"""
    Material structure
    Holds basic data about a material including name, density and
    composition in mass fraction.
"""
struct Material
    name::String
    density::Union{Missing, AbstractFloat}
    a::Dict{Int,AbstractFloat} # Optional: custom atomic weights for the keys in this Material
    massfraction::Dict{Int,AbstractFloat}
    Material(name::AbstractString, massfrac::Dict{Int,U}, density=missing, a::Dict{Int,V}=Dict{Int,Float64}()) where { U <: AbstractFloat, V <: AbstractFloat } =
        new(name, density, a, massfrac)
end

name(mat::Material) = mat.name

"""
    density(mat::Material)
    Density in g/cm^3
"""
density(mat::Material) = mat.density


"""
    material(name::AbstractString, massfrac::Dict{Element,<:AbstractFloat}, density=missing, a::Dict{Element,<:AbstractFloat}=Dict{Element,Float64}() )
    Construct a material from mass fractions and (optional) atomic weights.
"""
function material(name::AbstractString,
    massfrac::Dict{Element, U},
    density::Union{Missing,AbstractFloat}=missing,
    a::Dict{Element, V}=Dict{Element,Float64}()
) where { U <: AbstractFloat, V <: AbstractFloat }
    mf = Dict( (z(elm), val) for (elm, val) in massfrac )
    aw = Dict( (z(elm), val) for (elm, val) in a )
    return Material(name, mf, density, aw)
end

pure(elm::Element) =
    material("Pure "*symbol(elm), Dict{}(elm=>1.0), density(elm))

function Base.show(io::IO, mat::Material)
    print(io, name(mat)," = (")
    comma=""
    for (z, mf) in mat.massfraction
        print(io, @sprintf("%s%s = %0.4f", comma, element(z).symbol, mf))
        comma=", "
    end
    if !ismissing(mat.density)
        print(io, @sprintf(", %0.2f g/cc", mat.density))
    end
    print(io,")")
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

has(mat::Material, elm::Element) =
    haskey(mat.massfraction, z(elm))


"""
    atomicfraction(name::String, atomfracs::Dict{Element,Float64}, density = nothing, atomicweights::Dict{Element,Float64}=Dict())
Build a Material from atomic fractions (or stoichiometries).
"""
function atomicfraction(
    name::String,
    atomfracs::Dict{Element,U},
    density::Union{Missing,AbstractFloat} = missing,
    atomicweights::Dict{Element,V}=Dict{Element,Float64}() ) where { U <: Number, V <: AbstractFloat }
    aw(elm) = get(atomicweights, elm.number, a(elm))
    norm = sum(af*aw(elm) for (elm, af) in atomfracs)
    massfracs = Dict( (elm, af*aw(elm)/norm) for (elm, af) in atomfracs )
    return material(name, massfracs, density, atomicweights)
end

"""
    mac(mat::Material, energy::Float64)::Float64
Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, energy::Float64) =
    sum(mac(elm, energy) * massfraction(elm, mat) for elm in keys(mat))

"""
    mac(mat::Material, energy::Float64)::Float64
Compute the material MAC using the standard mass fraction weighted formula.
"""
mac(mat::Material, xray::CharXRay) =
    sum(mac(elm, xray) * massfraction(elm, mat) for elm in keys(mat))
