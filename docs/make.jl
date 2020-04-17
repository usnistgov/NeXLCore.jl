using Documenter
using NeXLCore
using Weave

rm(joinpath(@__DIR__,"build","gettingstarted.html"))
makedocs(modules = [NeXLCore], sitename = "NeXLCore")
weave(joinpath(@__DIR__,"..","examples","gettingstarted.jmd"), out_path=joinpath(@__DIR__,"build","gettingstarted.html"))

# deploydocs(repo = "github.com/JuliaLang/NeXLCore.jl.git")
