using .Gadfly

using Colors


const NeXLPalette = distinguishable_colors(
    66,
    Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)],
)[3:end]
const NeXLColorblind = distinguishable_colors(
    66,
    Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0), colorant"DodgerBlue4"],
    transform = deuteranopic,
)[3:end]

"""
    Gadfly.plot(transitions::AbstractVector{Transition}; mode=:Energy|:Weight, palette=NeXLPalette)

Plot either the :Energies or :Weights associated with the specified transitions over the range of supported elements.
"""
function Gadfly.plot(transitions::AbstractVector{Transition}; mode = :Energy, palette = NeXLPalette)
    if mode == :Energy
        plotXrayEnergies(transitions::AbstractVector{Transition}, palette = palette)
    elseif mode == :Weight
        plotXrayWeights(transitions::AbstractVector{Transition}, palette = palette)
    end
end

function plotXrayEnergies(transitions::AbstractVector{Transition}; palette = NeXLPalette)
    layers, names, colors = [], [], []
    for tr in transitions
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, tr)
                push!(x, z(elm))
                push!(y, energy(characteristic(elm, tr)))
            end
        end
        if !isempty(x)
            push!(names, repr(tr))
            push!(colors, NeXLPalette[length(colors)%length(palette)+1])
            push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[end])))
        end
    end
    println("Plotting $(length(layers)) layers.")
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Characteristic X-ray Energies"),
        # Gadfly.Guide.manual_color_key("Type", names, colors),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop),
    )
end

function plotXrayWeights(transitions::AbstractVector{Transition}; palette = NeXLPalette)
    layers, names, colors = [], [], []
    for tr in transitions
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, tr)
                push!(x, z(elm))
                push!(y, strength(elm, tr))
            end
            if !isempty(x)
                push!(names, repr(tr))
                push!(colors, palette[length(colors)%length(palette)+1])
                push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[end])))
            end
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Characteristic X-ray Energies"),
        Gadfly.Guide.manual_color_key("Type", name = names, colors = colors),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Line Weight"),
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop),
    )
end

"""
    Gadfly.plot(sss::AbstractVector{SubShell}; palette=NeXLPalette)

Plot the edge energies associated with the specified vector of SubShell objects.
"""
function Gadfly.plot(sss::AbstractVector{SubShell}; palette = NeXLPalette)
    plotEdgeEnergies(sss, palette = palette)
end


function plotEdgeEnergies(sss::AbstractVector{SubShell}; palette = NeXLPalette)
    layers, names, colors = [], [], []
    for sh in sss
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, sh)
                push!(x, z(elm))
                push!(y, energy(atomicsubshell(elm, sh)))
            end
            if !isempty(x)
                push!(names, repr(tr))
                push!(colors, palette[length(colors)%length(palette)+1])
                push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[end])))
            end
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Atomic Sub-Shell Energies"),
        Gadfly.Guide.manual_color_key("Type", name = names, colors = colors),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Edge Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop),
    )
end

"""
    compareMACs(elm::Element; palette=NeXLPalette)

Plot a comparison of the FFAST and Heinrich MAC tabulations for the specified Element.
"""
function compareMACs(elm::Element; palette = NeXLPalette)
    l1 = layer(ev -> log10(mac(elm, ev)), 100.0, 20.0e3, Geom.line, Gadfly.Theme(default_color = palette[1]))
    l2 = layer(ev -> log10(dtsamac(elm, ev)), 100.0, 20.0e3, Geom.line, Gadfly.Theme(default_color = palette[2]))
    Gadfly.plot(
        l1,
        l2,
        Gadfly.Guide.title("Comparing MACs"),
        Gadfly.Guide.manual_color_key("Type", ("Default", "Heinrich"), palette[1:2]),
        Gadfly.Guide.xlabel("Energy (eV)"),
        Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Coord.cartesian(xmin = 0.0, xmax = 20.0e3),
    )
end
