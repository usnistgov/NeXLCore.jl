using .MeshCat
using CoordinateTransformations

"""
    draw(vis::Visualizer, r::Region, colors::Dict{Material, Color}, name::String="Sample")

Draw the samples on the visualizer in translucent colors by Material. 
"""
function draw(vis::Visualizer, region::Region, colors::Dict{Material, Color}, name::String="Sample")
    for (i, ch) in enumerate(region.children)
        chname = "$name[$i]"
        setobject!(vis[chname], MeshObject(ch.shape, MeshBasicMaterial(color=RGBA(colors[ch.material],0.2f0))))
        settransform!(vis[chname],LinearMap([ 1.0e4 0.0 0.0; 0.0 1.0e4 0.0; 0.0 0.0 -1.0e4]))
        draw(vis, ch, colors, chname)
    end
end

function draw(vis::Visualizer, e0::Float64, sample::Region, num_trajectories::Int = 100)
    matcolors = colorize(sample)
    # Run many trajectories
    draw(vis, sample, matcolors, "Sample")
    for i in Base.OneTo(num_trajectories)
        pts, colors=Position[], RGB{Float32}[]
        trajectory(gun(Electron, e0, 1.0e-6), sample) do part, reg
            # Record the position and the color at each scatter point
            push!(pts, 10000.0*(position(part) .* (1.0e4, 1.0e4, -1.0e4)))
            push!(colors, matcolors[reg.material])
        end
        setobject!(vis["Trajectory[$i]"], PointCloud(pts, colors)) # Add replacements
    end
end