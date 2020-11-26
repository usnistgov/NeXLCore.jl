using LinearAlgebra
using StaticArrays

"""
`Position` : A point in 3-D.
"""
const Position = SVector{3,Float64}

"""
`Shape` is an abstract type that implements three methods.

    isinside(r::Shape, pos::Position)

Is `pos` strictly inside `r`?

    intersection(r::Shape, pos0::Particle, pos1::Particle)::Float64

Return a number `f` which represent the fraction of the distance between `pos0` and `pos1` that
first intersects the `Shape` `r`.  The intersection point will equal `pos0 .+ f*(pos1 .- pos0)`.
If `f` is between 0.0 and 1.0 then the intersection is on the interval between `pos0` and `pos1`.
If the ray from `pos0` towards `pos2` does not intersect `r` then this function returns Inf64.


    random_point_inside(sr::Shape)::Position

Generate a randomized point within the shape to permit a stochastic test of whether
a child region is fully contained within a parent region. 
"""
abstract type Shape end

"""
    RectangularShape <: Shape

A rectangular region oriented with the coordinate axes.
"""
struct RectangularShape <: Shape
    corner1::Position
    corner2::Position

    RectangularShape(c1::Position, c2::Position) = new(min.(c1,c2), max.(c1,c2))
end

Base.show(io::IO, rs::RectangularShape) = print(io, "Rect[$(rs.corner1), $(rs.corner2)]")

function isinside(rr::RectangularShape, pos::Position)
    all(pos[i] > rr.corner1[i] && pos[i] < rr.corner2[i] for i in eachindex(pos))
end

function intersection(rr::RectangularShape, pos1::Position, pos2::Position)::Float64
    _between(a, b, c)  = (a>b) && (a<c)
    t = Inf64
    for i in eachindex(pos1)
        j, k = i % 3 + 1, (i + 1) % 3 + 1
        if pos2[i] != pos1[i]
            u = (rr.corner1[i] - pos1[i]) / (pos2[i] - pos1[i])
            if (u > 0.0) && (u <= t) && #
                _between(pos1[j] + u * (pos2[j] - pos1[j]), rr.corner1[j], rr.corner2[j]) && # 
                _between(pos1[k] + u * (pos2[k] - pos1[k]), rr.corner1[k], rr.corner2[k])
                t = u
            end
            u = (rr.corner2[i] - pos1[i]) / (pos2[i] - pos1[i])
            if (u > 0.0) && (u <= t) && #
                _between(pos1[j] + u * (pos2[j] - pos1[j]), rr.corner1[j], rr.corner2[j]) && # 
                _between(pos1[k] + u * (pos2[k] - pos1[k]), rr.corner1[k], rr.corner2[k])
                t = u
            end
        end
    end
    return t
end

function random_point_inside(rr::RectangularShape)
    return Position( (rr.corner1 .+ rand(Float64, (3)) .* (rr.corner2 .- rr.corner1))...)
end


"""
    SphericalShape <: Shape

A spherical region.
"""
struct SphericalShape <: Shape
    center::Position
    radius::Float64

    SphericalShape(center::Position, radius::Float64) = new(center, abs(radius))
end

isinside(sr::SphericalShape, pos::Position) = norm(pos-sr.center) < sr.radius

function intersection(sr::SphericalShape, pos0::Position, pos1::Position)::Float64
    d, m = pos1 .- pos0, pos0 .- sr.center
    ma2, b = -2.0 * dot(d, d), 2.0 * dot(m, d)
    f = b^2 + ma2 * 2.0 * (dot(m, m) - sr.radius^2)
    if f >= 0.0 
       up, un = (b + sqrt(f)) / ma2, (b - sqrt(f)) / ma2
       return min(up < 0.0 ? Inf64 : up, un < 0.0 ? Inf64 : un)
    end
    return Inf64
end

function random_point_inside(sr::SphericalShape)
    p = 2.0 * rand(Float64, (3)) .- ( 1.0, 1.0 ,1.0 )
    while norm(p) >= 1.0
        p = 2.0 * rand(Float64, (3)) .- ( 1.0, 1.0 ,1.0 )
    end
    return Position( (sr.center + sr.radius*p)...)
end

"""
Particle represents a type that may be simulated using a transport Monte Carlo.  It must provide
these methods:

    position(el::Particle)::Position
    previous(el::Particle)::Position
    energy(el::Particle)::Float64

The position of the current and previous elastic scatter locations which are stored in that Particle type.

    T(prev::Position, curr::Position, energy::Energy) where {T <: Particle }
    T(el::T, 𝜆::Float64, 𝜃::Float64, 𝜑::Float64, ΔE::Float64) where {T <: Particle }

Two constructors: One to create a defined Particle and the other to create a new Particle based off
another which is translated by `λ` at a scattering angle (`θ`, `ϕ`) which energy change of `ΔE`

    transport(pc::T, mat::Material)::NTuple{4, Float64} where {T <: Particle }

A function that generates the values of ( `λ`, `θ`, `ϕ`, `ΔE`) for the specified `Particle` in the specified `Material`.
"""
abstract type Particle end

struct Electron <: Particle
    previous::Position
    current::Position
    energy::Float64 # eV

"""
    Electron(prev::Position, curr::Position, energy::Float64)
    Electron(el::Electron, 𝜆::Float64, 𝜃::Float64, 𝜑::Float64, ΔE::Float64)::Electron

Create a new `Electron` from this one in which the new `Electron` is a distance `𝜆` from the
first along a trajectory that is `𝜃` and `𝜑` off the current trajectory.
"""
    Electron(prev::Position, curr::Position, energy::Float64) = new(prev,curr,energy)   
    function Electron(el::Electron, 𝜆::Float64, 𝜃::Float64, 𝜑::Float64, ΔE::Float64)
        (u,v,w) = LinearAlgebra.normalize(position(el) .- previous(el))
        sc =  1.0 - abs(w) > 1.0e-8 ? #
            SA_F64[ #
                u*cos(𝜃)+sin(𝜃)*(u*w*cos(𝜑)-v*sin(𝜑))/sqrt(1.0 - w^2), #
                v*cos(𝜃)+sin(𝜃)*(u*w*cos(𝜑)+v*sin(𝜑))/sqrt(1.0 - w^2), #
                w*cos(𝜃)-sqrt(1.0 - w^2)*sin(𝜃)*cos(𝜑) # 
            ] : SA_F64[ #
                sign(w)*sin(𝜃)*cos(𝜑), #
                sign(w)*sin(𝜃)*sin(𝜑), #
                sign(w)*cos(𝜃) ]
        return new(el.current, el.current .+ 𝜆*sc, el.energy + ΔE)
    end
end

Base.show(io::IO, el::Electron) = print(io, "Electron[$current, $energy eV]")
position(el::Particle) = el.current
previous(el::Particle) = el.previous
energy(el::Particle) = el.energy

"""
    transport(pc::Electron, mat::Material, ecx=Liljequist1989, bethe=JoyLuo)::NTuple{4, Float64}

The default function defining elastic scattering and energy loss for an Electron.

Returns ( `λ`, `θ`, `ϕ`, `ΔE`) where `λ` is the mean path length, `θ` is the elastic scatter angle, `ϕ` is the azimuthal elastic scatter
angle and `ΔE` is the energy loss for transport over the distance `λ`.
"""
function transport(pc::Electron, mat::Material, ecx=Liljequist1989, bethe=JoyLuo)::NTuple{4, Float64}
    ( elm, 𝜆 ) = λ(ecx, mat, pc.energy)
    𝜆′=-𝜆*log(rand())
    𝜃 = acos(2.0*rand() - 1.0)
    dE = 𝜆′ * dEds(bethe, pc.energy, elm, mat[:Density]*mat[elm])
    return ( 𝜆′, rand(ecx, elm, pc.energy), 𝜃, dE ) 
end

"""
    gun(::Type{T}, energy::Float64, width::Float64=1.0e-7, initial::Position=SA_F64[0.0,0.0,-10.0], direction=SA_F64[0.0,0.0,1.0])::T where {T <: Particle}

A helper to construct the initial `Particle` in a randomized Gaussian distribution.
"""
function gun(::Type{T}, energy::Float64, width::Float64=1.0e-7, initial::Position=SA_F64[0.0,0.0,-10.0], direction=SA_F64[0.0,0.0,1.0])::T where {T <: Particle}
    r, th = sqrt(-2.0 * log(rand())) * width, 2.0π*rand()
    st = initial .+ SA_F64[ r*cos(th), r*sin(th), 0 ]
    return T(st .- direction, st, energy)
end

"""
    pathlength(el::Particle)

Length of the line segment represented by `el`.
"""
pathlength(el::Particle) = norm(position(el) .- previous(el))

intersection(r::Shape, p::Particle) = intersection(r, previous(p), position(p))


"""
    Region

A `Region` combines a `Shape` and a `Material` (with `:Density` property) and may fully contain zero or more child `Region`s.
"""
struct Region 
    shape::Shape
    material::Material
    parent::Union{Missing,Region}
    children::Vector{Region}

    function Region(sh::Shape, mat::Material, parent::Union{Missing,Region}, ntests=1000) 
        @assert mat[:Density] > 0.0
        res = new(sh, mat, parent, Region[])
        if !ismissing(parent)
            @assert all(_->isinside(parent.shape, random_point_inside(sh)), Base.OneTo(ntests)) "The child $sh is not fully contained within the parent $(parent.shape)."
            @assert all(ch->all(_->!isinside(ch.shape, random_point_inside(sh)), Base.OneTo(ntests)), parent.children) "The child $sh overlaps a child of the parent shape."
            push!(parent.children, res)
        end
        return res
    end
end

Base.show(io::IO, reg::Region) = print(io, "Region[$(reg.shape), $(reg.material), $(length(reg.children)) children]")

function parent(reg::Region)::Region
    ismissing(reg.parent) ? reg : reg.parent
end

"""
    childmost_region(reg::Region, pos::Position)

Find the inner-most `Region` within `reg` containing the point `pos`.
"""
function childmost_region(reg::Region, pos::Position)
    res = findfirst(ch->isinside(ch.shape,pos), reg.children)
    if !isnothing(res)
        return childmost_region(res, pos)
    end
    return reg
end

"""
    take_step(p::T, reg::Region, 𝜆::Float64, 𝜃::Float64, 𝜑::Float64)::Tuple{T, Region} where { T<: Particle}

Returns a `Tuple` containing a new `Particle` and the child-most `Region` in which the new `Particle` is found based
on a scatter event consisting a translation of up to `𝜆` mean-free path along a new direction given relative
to the current direction of `p` via the scatter angles `𝜃` and `𝜑`.

Returns the updated `Particle` reflecting the last trajectory step and the Region for the next step.
"""
function take_step(p::T, reg::Region, 𝜆::Float64, 𝜃::Float64, 𝜑::Float64, ΔE::Float64, ϵ::Float64 = 1.0e-12)::Tuple{T, Region} where {T <: Particle }
    newP, nextReg = T(p, 𝜆, 𝜃, 𝜑, ΔE), reg
    t = min( 
        intersection(reg.shape, newP), # Leave this Region
        (intersection(ch.shape, newP) for ch in reg.children)... # Enter a new child Region
    )
    if t <= 1.0
        newP = T(p, (t+ϵ)*𝜆, 𝜃, 𝜑, (t+ϵ)*ΔE)
        nextReg = childmost_region(parent(reg), position(newP))
    end
    return (newP, nextReg)
end

"""
    chamber(dim=100.0)

Construct a high vacuum chamber to act as the outer-most Region.
"""
chamber(dim=100.0) = #
    Region(RectangularShape(SA_F64[-dim,-dim,-dim],SA_F64[dim,dim,dim]), parse(Material,"H",density=5.0e-11), missing, 0)

"""
trajectory(p::T, reg::Region, scf::Function=transport, minE::Float64=50.0)

Run a single particle trajectory from `p` to `minE` or until the particle exits `reg`.

  * `eval(part::T, region::Region)` a function evaluated at each scattering point
  * `p` defines the initial position, direction and energy of the particle (often created with `gun(T, ...)`)
  * `reg` The outer-most region for the trajectory (usually created with `chamber()`)
  * `scf` A function from (<:Particle, Material) -> ( λ, θ, ϕ, ΔE ) that implements the transport dynamics
  * `minE` Stopping criterion
"""
function trajectory(eval::Function, p::T, reg::Region, scf::Function=(t::T, mat::Material) -> transport(t, mat), minE::Float64=50.0) where { T <: Particle }
    (pc, nextr) = (p, childmost_region(reg, current(p))) 
    while (pc.energy > minE) && isinside(reg.shape, position(pc))
        prevr = nextr
        (pc, nextr) = take_step(pc, nextr, scf(pc, nextr.material)...)
        eval(pc, prevr)
    end
end

