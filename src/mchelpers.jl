"""
    chamber(dim=100.0)

Construct a high vacuum chamber to act as the outer-most Region.
"""
chamber(dim = 100.0) = #
    Region(
        RectangularShape([-dim, -dim, -dim], 2.0 * [dim, dim, dim]),
        parse(Material, "H", density = 5.0e-11),
        nothing,
        "Chamber",
        0,
    )

"""
    gun(::Type{T}, energy::Float64, width::Float64=1.0e-7, initial::Position=Position(0.0,0.0,-10.0), direction=Position(0.0,0.0,1.0)::T where {T <: Particle}

A helper to construct the initial `Particle` in a randomized Gaussian distribution.
"""
function gun(
    ::Type{T},
    energy::Float64,
    width::Float64 = 1.0e-7,
    initial::Position = Position(0.0, 0.0, -10.0),
    direction = Position(0.0, 0.0, 1.0),
)::T where {T<:Particle}
    r, th = sqrt(-2.0 * log(rand())) * width, 2.0π * rand()
    st = initial .+ Position(r * cos(th), r * sin(th), 0.0)
    return T(st .- direction, st, energy)
end

"""
    bulk(mat::Material)

Construct a bulk homogeneous sample at the origin.
"""
function bulk(mat::Material)
    c = chamber()
    Region(RectangularShape([-0.1, -0.1, 0.0], [0.2, 0.2, 0.2]), mat, c, "Bulk homogeneous")
    return c
end

"""
    particle(mat::Material, radius::Float64, substrate::Union{Nothing,Material} = nothing)

Construct a spherical particle at the origin with an optional bulk substrate.
"""
function particle(
    mat::Material,
    radius::Float64,
    substrate::Union{Nothing,Material} = nothing,
)
    c = chamber()
    Region(SphericalShape([0.0, 0.0, -radius], radius), mat, c, "Particle")
    if !isnothing(substrate)
        Region(
            RectangularShape([-0.1, -0.1, 0.0], [0.2, 0.2, 0.2]),
            substrate,
            c,
            "Substrate",
        )
    end
    return c
end


"""
    coated_particle(mat::Material, radius::Float64, coating::Material, thickness::Float64, substrate::Union{Nothing,Material} = nothing)

Construct a coated particle on an optional substrate.

This model is a good example of how more complex models are constructed.  Use `dump(..)` to display the structure.  You'll notice that
the chamber serves as the root `Region`.  The `coating` and `substrate` are in the chamber (children of the chamber `Region`).  The
`particle` is a child of the `coating` `Region` because the particle is fully enclosed by the coating.  An electron inside of the coating will
enter the particle and leave the coating as soon as it enters the volume representing the particle.  The electron only appears to be
within the child-most region at any point in space.  So a typical trajectory might start in the chamber, enter the coating, traverse
the thickness of the coating and then enter the particle.  It may eventually leave the particle and reenter the coating, exit the coating
and reenter the chamber and finally come to rest in the substrate.
"""
function coated_particle(
    mat::Material,
    radius::Float64,
    coating::Material,
    thickness::Float64,
    substrate::Union{Nothing,Material} = nothing,
)
    c = chamber()
    cr = Region(
        SphericalShape([0.0, 0.0, -(radius + thickness)], radius + thickness),
        coating,
        c,
        "Coating",
    )
    Region(SphericalShape([0.0, 0.0, -(radius + thickness)], radius), mat, cr, "Particle")
    if !isnothing(substrate)
        Region(
            RectangularShape([-0.1, -0.1, 0.0], [0.2, 0.2, 0.2]),
            substrate,
            c,
            "Substrate",
        )
    end
    return c
end


"""
    thin_film(prs::Pair{Material, Float64}...; substrate::Union{Nothing,Material} = nothing)

Construct sample consisting of one or more thin films on an optional bulk substrate.
"""
function thin_film(
    prs::Pair{<:Material,Float64}...;
    substrate::Union{Nothing,Material} = nothing,
)
    c = chamber()
    t = 0.0
    for (i, (mat, thick)) in enumerate(prs)
        Region(RectangularShape([-0.1, -0.1, t], [0.2, 0.2, thick]), mat, c, "Layer[$i]")
        t += thick
    end
    if !isnothing(substrate)
        Region(
            RectangularShape([-0.1, -0.1, t], [0.2, 0.2, 0.2 - t]),
            substrate,
            c,
            "Substrate",
        )
    end
    return c
end

"""
    allmaterials(reg::Region)

Generates a list of all the unique materials in a sample Region. 
"""
function allmaterials(reg::Region)
    res = [reg.material]
    for ch in reg.children
        append!(res, allmaterials(ch))
    end
    return unique(res)
end

"""
    colorize(reg::Region)::Dict{<:Material, Color}

Generate a `Dict{Material, Color}` for all the `Materials` in the specified `Region`.
Designed for distinctive but not necessarily attractive colors.
"""
function colorize(reg::Region)::Dict{<:Material,Color}
    mats = allmaterials(reg)
    colors = distinguishable_colors(
        length(mats) + 2,
        [RGB(0.12, 0.14, 0.47), RGB(0.42, 0.63, 0.85), RGB(0.9, 0.9, 0.9)],
        transform = deuteranopic,
    )[3:end]
    return Dict(mat => col for (mat, col) in zip(mats, colors))
end

"""
    sample(ch::Region)

Extracts the sample portion from within the chamber.
"""
sample(ch::Region)::Vector{Region} = ch.children

function Base.dump(reg::Region, indent = 0)
    print("$(" → "^indent)$reg\n")
    for ch in reg.children
        dump(ch, indent + 1)
    end
end
