using .Gadfly

linecolors = ( "tomato", "darkseagreen4", "dodgerblue", "orange3", "goldenrod2",
    "gold4", "gold1", "olivedrab1", "deeppink4",  "cyan3", "chocolate",
    "darkcyan", "darksalmon", "blueviolet", "orchid4", "mediumpurple4",
    "maroon", "brown2", "chartreuse", "sienna", "firebrick" )

function plotXrayEnergies(transitions::AbstractArray{Transition})
    layers, names, colors = [], [], []
    for tr in transitions
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, tr)
                push!(x,z(elm))
                push!(y,energy(characteristic(elm,tr)))
            end
        end
        if !isempty(x)
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
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop))
end

function plotXrayWeights(transitions::AbstractArray{Transition})
    layers, names, colors = [], [], []
    for tr in transitions
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, tr)
                push!(x,z(elm))
                push!(y, strength(characteristic(elm,tr)))
            end
            if !isempty(x)
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
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop))
end

function plotEdgeEnergies(shells::AbstractArray{Shell})
    layers, names, colors = [], [], []
    for sh in shells
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, sh)
                push!(x,z(elm))
                push!(y,energy(atomicshell(elm,sh)))
            end
            if !isempty(x)
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
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop))
end


function compareMACs(elm::Element)
    l1 = layer(ev->log10(mac(elm,ev)), 100.0, 20.0e3, Geom.line, Gadfly.Theme(default_color="tomato"))
    l2 = layer(ev->log10(dtsamac(elm,ev)), 100.0, 20.0e3, Geom.line, Gadfly.Theme(default_color="darkseagreen4"))
    Gadfly.plot( l1, l2,
        Gadfly.Guide.title("Comparing MACs"),
        Gadfly.Guide.manual_color_key("Type", ("Default", "Heinrich"), ("tomato", "darkseagreen4")),
        Gadfly.Guide.xlabel("Energy (eV)"), Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Coord.cartesian(xmin = 0.0, xmax = 20.0e3))
end
