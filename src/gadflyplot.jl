using .Gadfly

linecolors = ( "tomato", "darkseagreen4", "dodgerblue", "orange3", "goldenrod2",
    "gold4", "gold1", "olivedrab1", "deeppink4",  "cyan3", "chocolate",
    "darkcyan", "darksalmon", "blueviolet", "orchid4", "mediumpurple4",
    "maroon", "brown2", "chartreuse", "sienna", "firebrick" )

function plotXrayEnergies(transitions::AbstractArray{Transition})
    layers, names, colors = [], [], []
    for tr in transitions
        x, y = [], []
        for elm in element.(80:92)
            if has(elm, tr)
                push!(x,z(elm))
                push!(y,energy(elm,tr))
            end
        end
        if length(x)>0
            push!(names, repr(tr))
            push!(colors, linecolors[length(colors) % length(linecolors) + 1])
            push!(layers, Gadfly.layer(x=x, y=y, Geom.point,
                Gadfly.Theme(default_color = colors[end] )))
        end
    end
    println("Plotting $(length(layers)) layers.")
    Gadfly.plot(layers...,
        Gadfly.Guide.title("Characteristic X-ray Energies"),
        # Gadfly.Guide.manual_color_key("Type", names, colors),
        Gadfly.Guide.xlabel("Atomic Number"), Guide.ylabel("Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = 80, xmax = 92))
end

function plotXrayWeights(transitions::AbstractArray{Transition})
    layers, names, colors = [], [], []
    for tr in transitions
        x, y = [], []
        for z in 1:92
            if has(element(z), tr)
                push!(x, z(element))
                push!(y, strength(element(z),tr))
            end
            if length(x)>0
                push!(names, repr(tr))
                push!(colors, linecolors[length(colors) % length(linecolors) + 1])
                push!(layers, Gadfly.layer(x=x, y=y, Geom.point,
                    Gadfly.Theme(default_color = colors[end] )))
            end
        end
    end
    Gadfly.plot(layers...,
        Gadfly.Guide.title("Characteristic X-ray Energies"),
        Gadfly.Guide.manual_color_key("Type", name = names, colors = colors),
        Gadfly.Guide.xlabel("Atomic Number"), Guide.ylabel("Line Weight"),
        Gadfly.Coord.cartesian(xmin = 80, xmax = 92))
end

function plotEdgeEnergies(shells::AbstractArray{Shell})
    layers, names, colors = [], [], []
    for sh in shells
        x, y = [], []
        for z in 1:92
            if has(element(z), sh)
                push!(x,z)
                push!(y,energy(element(z),sh))
            end
            if length(x)>0
                push!(names, repr(tr))
                push!(colors, linecolors[length(colors) % length(linecolors) + 1])
                push!(layers, Gadfly.layer(x=x, y=y, Geom.point,
                    Gadfly.Theme(default_color = colors[end] )))
            end
        end
    end
    Gadfly.plot(layers...,
        Gadfly.Guide.title("Atomic Shell Energies"),
        Gadfly.Guide.manual_color_key("Type", name = names, colors = colors),
        Gadfly.Guide.xlabel("Atomic Number"), Guide.ylabel("Edge Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = 1, xmax = 92))
end
