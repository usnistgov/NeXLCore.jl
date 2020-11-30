"""
    chamber(dim=100.0)

Construct a high vacuum chamber to act as the outer-most Region.
"""
chamber(dim=100.0) = #
    Region(RectangularShape([-dim,-dim,-dim],2.0*[dim,dim,dim]), parse(Material,"H",density=5.0e-11), missing, 0)

"""
    gun(::Type{T}, energy::Float64, width::Float64=1.0e-7, initial::Position=Position(0.0,0.0,-10.0), direction=Position(0.0,0.0,1.0)::T where {T <: Particle}

A helper to construct the initial `Particle` in a randomized Gaussian distribution.
"""
function gun(::Type{T}, energy::Float64, width::Float64=1.0e-7, initial::Position=Position(0.0,0.0,-10.0), direction=Position(0.0,0.0,1.0))::T where {T <: Particle}
    r, th = sqrt(-2.0 * log(rand())) * width, 2.0π*rand()
    st = initial .+ Position( r*cos(th), r*sin(th), 0.0 )
    return T(st .- direction, st, energy)
end

"""
    bulk(mat::Material)

Construct a bulk homogeneous sample at the origin.
"""
function bulk(mat::Material)
    c = chamber()
    Region(RectangularShape([-1.0e-2,-1.0e-2,0.0],[2.0e-2,2.0e-2,2.0e-2]), mat, c)
    return c
end

"""
    particle(mat::Material, radius::Float64, substrate::Union{Nothing,Material} = nothing)

Construct a spherical particle at the origin with an optional bulk substrate.
"""
function particle(mat::Material, radius::Float64, substrate::Union{Nothing,Material} = nothing)
    c = chamber()
    Region(SphericalShape([0.0, 0.0, -radius], radius), mat, c)
    if !isnothing(substrate)
        Region(RectangularShape([-1.0e-3, -1.e-30, 0.0],[2.0e-3,2.0e-3,2.0e-2]), substrate, c)
    end
    return c
end


"""
    coated_particle(mat::Material, radius::Float64, coating::Material, thickness::Float64, substrate::Union{Nothing,Material} = nothing)

Construct a coated particle on an optional substrate.
"""
function coated_particle(mat::Material, radius::Float64, coating::Material, thickness::Float64, substrate::Union{Nothing,Material} = nothing)
    c = chamber()
    cr = Region(SphericalShape([0.0, 0.0, -(radius+thickness)], radius+thickness), coating, c)
    Region(SphericalShape([0.0, 0.0, -(radius+thickness)], radius), mat, cr)
    if !isnothing(substrate)
        Region(RectangularShape([-1.0e-3, -1.0e-3, 0.0],[2.0e-3,2.0e-3,2.0e-2]), substrate, c)
    end
    return c
end


"""
    thinfilm(prs::Pair{Material, Float64}...; substrate::Union{Nothing,Material} = nothing)

Construct sample consisting of one or more thin films on an optional bulk substrate.
"""
function thinfilm(prs::Pair{Material, Float64}...; substrate::Union{Nothing,Material} = nothing)
    c = chamber()
    t = 0.0
    for (mat, thick) in prs
        Region(RectangularShape([-1.0e-3,-1.0e-3,t],[2.0e-3,2.0e-3,t+thick]), mat, c)
        t+=thick
    end
    if !isnothing(substrate)
        Region(RectangularShape([-1.0e-3, -1.0e-3, t],[2.0e-3,2.0e-3,2.0e-2-t]), substrate, c)
    end
    return c
end

"""
    allmaterials(reg::Region)

Generates a list of all the unique materials in a sample Region. 
"""
function allmaterials(reg::Region)
    res = [ reg.material ]
    for ch in reg.children
        append!(res, allmaterials(ch))
    end
    return unique(res)
end

"""
    colorize(reg::Region)::Dict{Material, Color}

Generate a `Dict{Material, Color}` for all the `Materials` in the specified `Region`.
Designed for distinctive but not necessarily attractive colors.
"""
function colorize(reg::Region)::Dict{Material, Color}
    mats = allmaterials(reg)
    colors = distinguishable_colors(
        length(mats)+2,
        [RGB(0.12, 0.14, 0.47), RGB(0.42, 0.63, 0.85), RGB(0.9,0.9,0.9)],
        transform = deuteranopic,
    )[3:end]
    return Dict( mat=>col for (mat,col) in zip(mats, colors))
end

"""
    sample(ch::Region)

Extracts the sample portion from within the chamber.
"""
sample(ch::Region)::Vector{Region} = ch.children
