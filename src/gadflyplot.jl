using .Gadfly

using Colors
using Pkg.Artifacts
using CSV
using DataFrames
using Statistics

const NeXLPalette = distinguishable_colors(
    66,
    Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)],
    transform = deuteranopic,
)[3:end]
const NeXLColorblind = distinguishable_colors(
    66,
    Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0), colorant"DodgerBlue4"],
    transform = deuteranopic,
)[3:end]

#Compose.parse_colorant(c::Array{<:Colorant,1}) = c

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
    layers, names = [], String[]
    colors = distinguishable_colors(length(transitions) + 2, Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)])
    for (i, tr) in enumerate(transitions)
        x, y = [], []
        for elm in element.(elementrange())
            if has(elm, tr)
                push!(x, z(elm))
                push!(y, energy(characteristic(elm, tr)))
            end
        end
        if !isempty(x)
            push!(names, repr(tr))
            append!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Characteristic X-ray Energies"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = elementrange().start, xmax = elementrange().stop),
    )
end

function plotXrayWeights(transitions::AbstractVector{Transition}, schoonjan::Bool = false)
    layers, names = [], String[]
    colors = distinguishable_colors(length(transitions) + 2, Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)])
    for (i, tr) in enumerate(transitions)
        x, y = [], []
        for elm in element.(elementrange())
            if has(elm, tr)
                push!(x, z(elm))
                push!(y, schoonjan ? normweight(characteristic(elm, tr)) : strength(characteristic(elm, tr)))
            end
        end
        if !isempty(x)
            push!(names, repr(tr))
            append!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    if schoonjan  # Compare to Shoonjan's xraylib on GitHub
        function parse2(tr)
            ss1, ss2 = missing, missing
            try
                if tr[1] == 'K'
                    ss1 = n"K1"
                    ss2 = (tr[2] ≠ 'O' && tr[2] ≠ 'P') ? parse(SubShell, tr[2:end]) : (tr[2] == 'O' ? n"O3" : n"P3")
                else
                    ss1 = parse(SubShell, tr[1:2])
                    ss2 = (tr[3] ≠ 'O' && tr[3] ≠ 'P') ? parse(SubShell, tr[3:end]) : (tr[3] == 'O' ? n"O3" : n"P3")
                end
                return exists(ss1, ss2) ? Transition(ss1, ss2) : missing
            catch
                return missing
            end
        end
        artpath = downloadschoonjan()
        csv = CSV.read("$artpath\\radrate.dat", delim = ' ', ignorerepeated = true, header = 0) |> DataFrame
        insertcols!(csv, 3, Column2p = parse2.(csv[!, :Column2]))
        filter!(r -> (!ismissing(r[:Column2p])) && (r[:Column2p] in transitions), csv)
        insertcols!(csv, 3, Column2pp = CategoricalArray(map(n -> "Schoonj $n", csv[!, :Column2p])))
        append!(layers, layer(csv, x = :Column1, y = :Column3, color = :Column2pp))
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Characteristic X-ray Weights"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Weight"),
        Gadfly.Coord.cartesian(xmin = elementrange().start, xmax = elementrange().stop),
    )
end

"""
    Gadfly.plot(sss::AbstractVector{SubShell}, mode=:EdgeEnergy|:FluorescenceYield; palette=NeXLPalette)

Plot the edge energies/fluorescence yields associated with the specified vector of SubShell objects.
"""
function Gadfly.plot(sss::AbstractVector{SubShell}, mode = :EdgeEnergy)
    if mode == :FluorescenceYield
        plotFluorescenceYield(sss::AbstractVector{SubShell})
    else
        plotEdgeEnergies(sss)
    end
end

# Download the data once per installation...
function downloadschoonjan()
    artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")
    hash = artifact_hash("schoonjan", artifacts_toml)
    if hash == nothing || !artifact_exists(hash)
        hash = create_artifact() do artifact_dir
            # We create the artifact by simply downloading a few files into the new artifact directory
            download("https://github.com/tschoonj/xraylib/blob/master/data/fluor_yield.dat?raw=true", joinpath(artifact_dir, "fluor_yield.dat"))
            download("https://github.com/tschoonj/xraylib/blob/master/data/radrate.dat?raw=true", joinpath(artifact_dir, "radrate.dat"))
        end
        bind_artifact!(artifacts_toml, "schoonjan", hash)
        @info "Downloaded Schoonjan data into Archive."
    end
    return artifact_path(hash)
end

function plotFluorescenceYield(sss::AbstractVector{SubShell}, schoonjan::Bool = false)
    layers, names = [], String[]
    colors = distinguishable_colors(length(sss) + 2, Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)])
    for (i, sh) in enumerate(sss)
        x, y = [], []
        for elm in element.(elementrange())
            if has(elm, sh)
                push!(x, z(elm))
                push!(y, fluorescenceyield(atomicsubshell(elm, sh), NeXL))
            end
        end
        if !isempty(x)
            push!(names, repr(sh))
            append!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    if schoonjan
        artpath = downloadschoonjan()
        csv = CSV.read("$artpath\\fluor_yield.dat", delim = ' ', ignorerepeated = true, header = 0) |> DataFrame
        sssname = [repr(ss) for ss in sss]
        filter!(r -> r[:Column2] in sssname, csv)
        insertcols!(csv, 3, Column2p = CategoricalArray(map(n -> "Schoonj $n", csv[!, :Column2])))
        append!(layers, layer(csv, x = :Column1, y = :Column3, color = :Column2p))
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Fluourescence Yield"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Yield (Fractional)"),
        Scale.y_log10(maxvalue = 1.0),
        Gadfly.Coord.cartesian(xmin = elementrange().start, xmax = elementrange().stop),
    )
end


function plotEdgeEnergies(sss::AbstractVector{SubShell})
    layers, names = [], String[]
    colors = distinguishable_colors(length(sss) + 2, Color[RGB(253 / 255, 253 / 255, 241 / 255), RGB(0, 0, 0)])
    for (i, sh) in enumerate(sss)
        x, y = [], []
        for elm in element.(elementrange())
            if has(elm, sh)
                push!(x, z(elm))
                push!(y, energy(atomicsubshell(elm, sh)))
            end
        end
        if !isempty(x)
            push!(names, repr(sh))
            append!(layers, Gadfly.layer(x = x, y = y, Geom.point, Gadfly.Theme(default_color = colors[i+2])))
        end
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.title("Atomic Sub-Shell Energies"),
        Gadfly.Guide.manual_color_key("Type", names, colors[3:end]),
        Gadfly.Guide.xlabel("Atomic Number"),
        Guide.ylabel("Edge Energy (eV)"),
        Gadfly.Coord.cartesian(xmin = elementrange().start, xmax = elementrange().stop),
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
        Gadfly.Guide.manual_color_key("Type", [ "Default/FFAST", "Heinrich" ], palette[1:2]),
        Gadfly.Guide.xlabel("Energy (eV)"),
        Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Coord.cartesian(xmin = 0.0, xmax = 20.0e3),
    )
end

"""
    plot(alg::Type{<:NeXLAlgorithm}, elm::Union{Element,Material}; palette = NeXLPalette, xmax=20.0e3)

Plot a MAC tabulations for the specified Element or Material.
"""
function Gadfly.plot(alg::Type{<:NeXLAlgorithm}, elm::Union{Element,Material}; palette = NeXLPalette, xmax = 20.0e3)
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
function Gadfly.plot(
    alg::Type{<:NeXLAlgorithm},
    elms::AbstractVector;
    palette = NeXLPalette,
    xmin = 100.0,
    xmax = 20.0e3,
)
    layers, colors, names = Layer[], Color[], String[]
    for (i, elm) in enumerate(elms)
        append!(
            layers,
            layer(ev -> log10(mac(elm, ev, alg)), xmin, xmax, Geom.line, Gadfly.Theme(default_color = palette[i])),
        )
        push!(colors, palette[i])
        push!(names, name(elm))
    end
    Gadfly.plot(
        layers...,
        Gadfly.Guide.xlabel("Energy (eV)"),
        Guide.ylabel("log₁₀(MAC (cm²/g))"),
        Gadfly.Guide.manual_color_key("Material", names, colors),
        Gadfly.Coord.cartesian(xmin = max(0.0, xmin - 100.0), xmax = xmax),
    )
end

function Gadfly.plot(mats::AbstractVector{Material}; known::Union{Material, Missing}=missing, delta::Bool=false, label::AbstractString="Material")
    allelms = collect(union(map(keys,mats)...))
    xs = [ name(mat) for mat in mats ]
    layers, names, colors = Layer[], String[], RGB{Float64}[]
    if ismissing(known)
        known = material("the Mean", Dict{Element,Float64}( elm=>mean([value(mat[elm]) for mat in mats]) for elm in allelms))
    end
    for (i, elm) in enumerate(allelms)
        if delta
            append!(layers,
                layer(x=xs, y=[ value(mat[elm])-known[elm] for mat in mats],
                    ymin = [ value(mat[elm])-σ(mat[elm])-known[elm] for mat in mats],
                    ymax = [ value(mat[elm])+σ(mat[elm])-known[elm] for mat in mats],
                    Gadfly.Theme(default_color = NeXLPalette[i]), Geom.errorbar, Geom.point
                )
            )
        else
            append!(layers,
                layer(x=xs, y=[ value(mat[elm]) for mat in mats],
                    ymin = [ value(mat[elm])-σ(mat[elm]) for mat in mats],
                    ymax = [ value(mat[elm])+σ(mat[elm]) for mat in mats],
                    Gadfly.Theme(default_color = NeXLPalette[i]), Geom.errorbar, Geom.point
                )
            )
        end
        push!(names, name(elm))
        push!(colors, NeXLPalette[i])
    end
    lighten(col)=weighted_color_mean(0.2, RGB(col), colorant"white")
    if delta
        plot(layers..., Guide.ylabel("Δ(Mass Fraction)"),
            Guide.xlabel(label),
            Guide.manual_color_key("Element", names, colors),
            Guide.title("Difference from $(known)"),
            Geom.hline(color="black"), yintercept=[0.0])
    else
        plot(layers..., Guide.ylabel("Mass Fraction"), Guide.xlabel(label), #
                Guide.manual_color_key("Element", names, colors),
                yintercept=[known[elm] for elm in allelms],
                Geom.hline(color=[ lighten(col) for col in colors], style=:dash))
    end
end

function plot2(mats::AbstractVector{Material}; known::Union{Material, Missing}=missing, label::AbstractString="Material")
    allelms = sort(convert(Vector{Element},collect(union(map(keys,mats)...))))
    elmcol = Dict(elm=>NeXLPalette[i] for (i, elm) in enumerate(allelms))
    xs, ymin, ymax, ygroups, colors = String[], Float64[], Float64[], Element[], Color[]
    for mat in mats
        append!(xs, [ name(mat) for elm in keys(mat) ])
        append!(ymin, [ value(mat[elm])-σ(mat[elm]) for elm in keys(mat) ])
        append!(ymax, [ value(mat[elm])+σ(mat[elm]) for elm in keys(mat) ])
        append!(colors, [elmcol[elm] for elm in keys(mat)])
        append!(ygroups, collect(keys(mat)))
    end
    plot(x=xs, ymin=ymin, ymax=ymax, color=colors, ygroup=ygroups,
        Geom.subplot_grid(Geom.errorbar, free_y_axis=true), Scale.ygroup(labels=elm->symbol(elm), levels=allelms),
        Guide.xlabel(label), Guide.ylabel("Mass Fraction by Element"))
end
