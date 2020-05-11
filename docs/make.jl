using Documenter
using NeXLCore
using Weave

rm(joinpath(@__DIR__,"build","gettingstarted.html"))
makedocs(modules = [NeXLCore], sitename = "NeXLCore")
# deploydocs(repo = "github.com/JuliaLang/NeXLCore.jl.git")
