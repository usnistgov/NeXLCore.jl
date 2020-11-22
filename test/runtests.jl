using Test
using Random
using NeXLCore

@testset "NeXLCore Tests" begin
    include("xray.jl")
    include("material.jl")
    include("stoichiometry.jl")
    include("matu.jl")
    include("kratio.jl")
    include("eta.jl")
    include("custommacs.jl")
    include("bethe.jl")
    include("meanionizationpotential.jl")
    include("electron.jl")
    include("scatter.jl")
end

# To test coverage (from NWMR's Win10 box)
# 1. Open the command line
# 2. cd c:\Users\nritchie\.julia\dev\NeXLCore\test
# 3. > C:\Users\nritchie\AppData\Local\Programs\Julia\Julia-1.4.1\bin\julia.exe --code-coverage=user
# 4. julia> using Pkg
# 5. julia> Pkg.test("NeXLCore",coverage=true)
# 6. Look at .cov files in src for test coverage data
