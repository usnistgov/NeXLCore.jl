
"""
Stage coordinates are held in a `Dict{Symbol,Float64}`.  The nominal units
are mm and radians.  The default labels are `:X`, `:Y`, `:Z`, `:R`, `:T`, and `:B`.

The `StageMapping` abstract class defines the relationship between stage coordinates
and image (or hyperimage) pixels through the functions `image2stage(...)` and
`stage2image(...)`.

Since every instrument is different, you will likely have to implement your own 
`StageMapping` struct with its own functions.  `DefaultStageMapping` is provided
as an example of a simple mapping.
"""

abstract type StageMapping end

"""
The `DefaultStageMapping` assumes a matching Cartesian coordinate system for both the stage and image.  It only assumes
`:X`, `:Y` stage motion.
"""
struct DefaultStageMapping <: StageMapping end

"""
    image2stage(::Type{DefaultStageMapping}, stage_coord::Dict{Symbol,Float64}, img_coord::Dict{Symbol,Float64}, theta::Float64)

Given the stage coordinate at the center of the image `stage_coord`, the pixel coordinate of interest (in the same units as
`stage_coord`) and the image rotation `theta`, compute the stage coordinate that would bring the pixel into the center of the
image. 
"""
function image2stage(::Type{DefaultStageMapping}, stage_coord::Dict{Symbol,Float64}, img_coord::Dict{Symbol,Float64}, theta::Float64)
    rot = [ cos(theta) -sin(theta) ; sin(theta) cos(theta) ]
    xy = [ stage_coord[:X], stage_coord[:Y] ] + rot * [ img_coord[:X], img_coord[:Y] ]
    return Dict{Symbol,Float64}(:X=>xy[1], :Y=>xy[2])
end

"""
    stage2image(::Type{DefaultStageMapping}, stage_coord::Dict{Symbol,Float64}, centered_coord::Dict{Symbol,Float64}, theta::Float64)

`stage2image(...)` is the inverse function of `image2stage(...)`.  Given the stage coordinate of the center of the image 
`stage_coord` and the stage coordinate that would center the pixel of interest `centered_coord`, compute the pixel coordinate
corresponding to the `centered_coordinate` when the stage is at `stage_coord`.
"""
function stage2image(::Type{DefaultStageMapping}, stage_coord::Dict{Symbol,Float64}, centered_coord::Dict{Symbol,Float64}, theta::Float64)
    roti = [ cos(theta) sin(theta) ; -sin(theta) cos(theta) ]
    xy = roti * ([ centered_coord[:X], centered_coord[:Y] ] - [ stage_coord[:X], stage_coord[:Y] ])
    return Dict{Symbol,Float64}(:X=>xy[1], :Y=>xy[2])
end
