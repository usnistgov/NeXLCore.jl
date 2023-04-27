using Test
using NeXLCore

@testset "NeXLCore Tests" begin
    include("bethe.jl")
    include("compton.jl")
    include("custommacs.jl")
    include("electron.jl")
    include("element.jl")
    include("eta.jl")
    include("kratio.jl")
    include("mac.jl")
    include("material.jl")
    include("matparser.jl")
    include("stoichiometry.jl")
    include("matu.jl")
    include("meanionizationpotential.jl")
    include("scattering.jl")
    include("staging.jl")
    include("standardize.jl")
    include("stoichiometry.jl")
    include("xray.jl")
    include("mats.jl")
end

# To test coverage (from NWMR's Win10 box)
# 1. Open the command line
# 2. cd c:\Users\nritchie\.julia\dev\NeXLCore\test
# 3. > C:\Users\nritchie\AppData\Local\Programs\Julia\Julia-1.4.1\bin\julia.exe --code-coverage=user
# 4. julia> using Pkg
# 5. julia> Pkg.test("NeXLCore",coverage=true)
# 6. Look at .cov files in src for test coverage data
