using Documenter
using NeXLCore


function addNISTHeaders(htmlfile::String)
    # read HTML
    html = transcode(String,read(htmlfile))
    # Find </head>
    i = findfirst(r"</[Hh][Ee][Aa][Dd]>", html)
    # Already added???
    j = findfirst("nist-header-footer", html)
    if isnothing(j) && (!isnothing(i))
        # Insert nist-pages links right before </head>
        res = html[1:i.start-1]*
            "<link rel=\"stylesheet\" href=\"https://pages.nist.gov/nist-header-footer/css/nist-combined.css\">\n"*
            "<script src=\"https://pages.nist.gov/nist-header-footer/js/jquery-1.9.0.min.js\" type=\"text/javascript\" defer=\"defer\"></script>\n"*
            "<script src=\"https://pages.nist.gov/nist-header-footer/js/nist-header-footer.js\" type=\"text/javascript\" defer=\"defer\"></script>\n"*
            html[i.start:end]
        write(htmlfile, res)
        println("Inserting NIST header/footer into $htmlfile")
    end
    return htmlfile
end

include(joinpath("..","weave","buildweave.jl"))

pages = [
            # "Start Page" => "index.md",
            "Getting Started" => "gettingstarted.md",
            "Materials" => "material.md",
            "Bremsstrahlung" => "bremsstrahlung.md",
            "Mass Fraction Uncertainty" => "Au60Ag40unc.md",
            "Backscatter" => "eta.md",
            "Fluorescence Yield" => "fluoryield.md",
            "O by Stoichiometry" => "OByStoic.md",
            "Mean Ionization Potential" => "meanionizationpotential.md",
            "API: Structures and Methods" => "methods.md"
         ]

makedocs(
    doctest=false,
    sitename = "NeXLCore",
    modules = [ NeXLCore ],
    pages = pages
)

addNISTHeaders(joinpath(@__DIR__, "build","index.html"))
addNISTHeaders.(map(name->joinpath(@__DIR__, "build", splitext(name)[1], "index.html"), map(p->p.second, pages)))

rm.(joinpath(@__DIR__,"src",page[2]) for page in pages[1:end-1])

# deploydocs(repo = "github.com/JuliaLang/NeXLCore.jl.git")
