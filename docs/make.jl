using Documenter
using NeXLCore

include(joinpath("..","weave","buildweave.jl"))

pages = [
            "Getting Started" => "gettingstarted.md",
            "Bremsstrahlung" => "bremsstrahlung.md",
            "Mass Fraction Uncertainty" => "Au60Ag40unc.md",
            "Backscatter" => "eta.md",
            "Fluorescence Yield" => "fluoryield.md",
            "O by Stoichiometry" => "OByStoic.md",
            "Mean Ionization Potential" => "meanionizationpotential.md"
         ]

makedocs(
    doctest=false,
    sitename = "NeXLCore",
    modules = [ NeXLCore ],
    pages = pages
)

rm.(joinpath("src",page[2]) for page in pages)

# deploydocs(repo = "github.com/JuliaLang/NeXLCore.jl.git")
