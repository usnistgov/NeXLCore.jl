using Weave

let start_dir = pwd()
    cd(@__DIR__)
    outpath = normpath(joinpath(@__DIR__, "..", "docs", "src"))
    @show outpath
    if !isdirpath(outpath)
        mkpath(outpath)
    end

    weaveit(name) = weave(
        joinpath(@__DIR__, "$name.jmd"),
        out_path = joinpath(outpath, "$name.md"),
        doctype = "github",
    )

    weaveit.((
        "gettingstarted",
        "OByStoic",
        "fluoryield",
        "meanionizationpotential",
        "Au60Ag40unc",
        "eta",
        "bremsstrahlung",
        "material"
    ))

    cd(start_dir)
end
