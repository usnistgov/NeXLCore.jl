using Test
using NeXLCore
using Statistics

@testset "Materials Struct" begin

    mats = Materials("TEST", [ n"Ar", n"Tc", n"Pu", n"Ac" ], Float32, (100, 100); atomicweights = Dict(n"Tc"=>99.0 ), properties = Dict{Symbol,Any}(:Density=>2.2))

    ms = reshape( Base.similar(mat"(0.3±0.01)*Ar+(0.2±0.01)*Tc+(0.4±0.01)*Pu+(0.1±0.01)*Ac", 100*100), ( 100, 100 ))

    foreach(ci->mats[ci] = ms[ci], CartesianIndices(mats))

    @test isapprox(mean(mats[n"Tc"]), 0.2f0, atol=0.002f0)   
    @test isapprox(mean(mats[n"Ar"]), 0.3f0, atol=0.002f0)    
    @test isapprox(mean(mats[n"Pu"]), 0.4f0, atol=0.002f0)

    @test size(mats) == (100, 100)
    @test ndims(mats) == 2
    @test a(n"Pu", mats)==244.0
    @test a(n"Tc", mats)==99.0
    @test properties(mats)[:Density]==2.2
    @test length(properties(mats))==1

    nmats=asnormalized(mats)
    @test isapprox(mean(analyticaltotal.(nmats)), 1.0, atol=1.0e-8)
    @test n"Pu" in keys(nmats)
    @test !(n"Fe" in keys(mats))
    @test keys(mats) == keys(nmats)
end