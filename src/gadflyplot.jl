using .Gadfly

using Colors

const NeXLPalette = ( # This palette - https://flatuicolors.com/palette/nl
	RGB(234/255, 32/255, 39/255), RGB(27/255, 20/255, 100/255),
	RGB(0/255, 98/255, 102/255), RGB(87/255, 88/255, 187/255),
	RGB(11/2551, 30/255, 81/255), RGB(247/255, 159/255, 31/255),
	RGB(18/255, 137/255, 167/255), RGB(163/255, 203/255, 56/255),
	RGB(217/255, 128/255, 250/255), RGB(181/255, 52/255, 113/255),
	RGB(238/255, 90/255, 36/255), RGB(6/255, 82/255, 221/255),
	RGB(0/255, 148/255, 50/255), RGB(153/255, 128/255, 250/255),
	RGB(131/255, 52/255, 113/255), RGB(255/255, 195/255, 18/255),
	RGB(18/255, 203/255, 196/255), RGB(196/255, 229/255, 56/255),
	RGB(253/255, 167/255, 223/255), RGB(237/255, 76/255, 103/255) )

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
            push!(colors, NeXLPalette[length(colors) % length(NeXLPalette) + 1])
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
                push!(colors, NeXLPalette[length(colors) % length(NeXLPalette) + 1])
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
                push!(colors, NeXLPalette[length(colors) % length(NeXLPalette) + 1])
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
