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
        plotXrayEnergies(transitions)
    elseif mode == :Weight
        plotXrayWeights(transitions)
    end
end

function plotXrayEnergies(transitions::AbstractVector{Transition}; palette = NeXLPalette)
    layers, names = [], []
    colors = distinguishable_colors(
        length(transitions)+2,
        Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)],
    )
    for (i,tr) in enumerate(transitions)
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, tr)
                push!(x, z(elm))
                push!(y, energy(characteristic(elm, tr)))
            end
        end
        if !isempty(x)
            push!(names, repr(tr))
            push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Characteristic X-ray Energies"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop),
    )
end

function plotXrayWeights(transitions::AbstractVector{Transition})
    layers, names = [], []
    colors = distinguishable_colors(
        length(transitions)+2,
        Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0) ],
    )
    for (i,tr) in enumerate(transitions)
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, tr)
                push!(x, z(elm))
                push!(y, strength(characteristic(elm, tr)))
            end
        end
        if !isempty(x)
            push!(names, repr(tr))
            push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Characteristic X-ray Weights"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Weight"),
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop),
    )
end

"""
    Gadfly.plot(sss::AbstractVector{SubShell}, mode=:EdgeEnergy|:FluorescenceYield; palette=NeXLPalette)

Plot the edge energies/fluorescence yields associated with the specified vector of SubShell objects.
"""
function Gadfly.plot(sss::AbstractVector{SubShell}, mode=:EdgeEnergy)
    if mode==:FluorescenceYield
        plotFluorescenceYield(sss::AbstractVector{SubShell})
    else
        plotEdgeEnergies(sss)
    end
end

function plotFluorescenceYield(sss::AbstractVector{SubShell})
    layers, names = [], []
    colors = distinguishable_colors(
        length(sss)+2,
        Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)],
    )
    for (i, sh) in enumerate(sss)
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, sh)
                push!(x, z(elm))
                push!(y, fluorescenceyield(atomicsubshell(elm, sh)))
            end
        end
        if !isempty(x)
            push!(names, repr(sh))
            push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Fluourescence Yield"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Yield (Fractional)"),
        Scale.y_log10(maxvalue=1.0),
        Gadfly.Coord.cartesian(xmin = elementRange().start, xmax = elementRange().stop),
    )
end


function plotEdgeEnergies(sss::AbstractVector{SubShell})
    layers, names = [], []
    colors = distinguishable_colors(
        length(sss)+2,
        Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)],
    )
    for (i, sh) in enumerate(sss)
        x, y = [], []
        for elm in element.(elementRange())
            if has(elm, sh)
                push!(x, z(elm))
                push!(y, energy(atomicsubshell(elm, sh)))
            end
        end
        if !isempty(x)
            push!(names, repr(sh))
            push!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Atomic Sub-Shell Energies"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
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
    l1 = layer(ev -> log10(mac(elm, ev, FFASTDB)), 100.0, 20.0e3, Geom.line, Gadfly.Theme(default_color = palette[1]))
    l2 = layer(ev -> log10(mac(elm, ev, DTSA)), 100.0, 20.0e3, Geom.line, Gadfly.Theme(default_color = palette[2]))
    Gadfly.plot(
        l1,
        l2,
        Gadfly.Guide.title("MAC - $elm"),
        Gadfly.Guide.manual_color_key("Type", ("Default/FFAST", "Heinrich"), palette[1:2]),
        Gadfly.Guide.xlabel("Energy (eV)"),
        Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Coord.cartesian(xmin = 0.0, xmax = 20.0e3),
    )
end

"""
    plot(alg::Type{<:NeXLAlgorithm}, elm::Union{Element,Material}; palette = NeXLPalette, xmax=20.0e3)

Plot a MAC tabulations for the specified Element or Material.
"""
function Gadfly.plot(alg::Type{<:NeXLAlgorithm}, elm::Union{Element,Material}; palette = NeXLPalette, xmax=20.0e3)
    l1 = layer(ev -> log10(mac(elm, ev, alg)), 100.0, xmax, Geom.line, Gadfly.Theme(default_color = palette[1]))
    Gadfly.plot(
        l1,
        Gadfly.Guide.title("MAC - $(name(elm))"),
        Gadfly.Guide.xlabel("Energy (eV)"),
        Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Coord.cartesian(xmin = 0.0, xmax = xmax),
    )
end

"""
    compareMACs(elm::Element; palette=NeXLPalette)

Plot a comparison of the FFAST and Heinrich MAC tabulations for the specified Element or Material.
"""
function Gadfly.plot(alg::Type{<:NeXLAlgorithm}, elms::AbstractVector; palette = NeXLPalette, xmin=100.0, xmax=20.0e3)
    layers, colors, names = Layer[], Color[], String[]
    for (i,elm) in enumerate(elms)
        append!(layers, layer(ev -> log10(mac(elm, ev, alg)), xmin, xmax, Geom.line, Gadfly.Theme(default_color = palette[i])))
        push!(colors, palette[i])
        push!(names, name(elm))
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.xlabel("Energy (eV)"),
        Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Guide.manual_color_key("Material", names, colors),
        Gadfly.Coord.cartesian(xmin = max(0.0,xmin-100.0), xmax = xmax),
    )
end
