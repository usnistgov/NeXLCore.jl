using Test
using NeXLCore

@testset "Compton" begin
    cxr = n"U L3-M5"
    @test isapprox(comptonShift(π/4, cxr), 1.0/(1.0+(energy(cxr)/0.511e6)*(1.0-cos(π/4))), rtol=1.0e-5)
    @test isapprox(energy(cxr)*(1.0 - comptonShift(π, cxr)), 689.0, atol=1.0)
    @test isapprox(comptonEnergy(π/3, cxr), energy(cxr)/(1.0+(energy(cxr)/0.511e6)*(1.0-cos(π/3))), rtol=1.0e-5)
    @test isapprox(energy(cxr) - comptonEnergy(π, cxr), 689.0, atol=1.0)
    @test energy(cxr) - comptonEnergy(π/2, cxr) < 689.0
    # These test values are taken from Wikipedia
    @test isapprox(comptonAngular(π/2, 2.75) / comptonAngular(0.0, 2.75), 0.5, atol=0.0001)
    @test isapprox(comptonAngular(π/2, 60.0e3) / comptonAngular(0.0, 60.0e3), 0.4, atol=0.02)
    @test isapprox(comptonAngular(π/2, 511.0e3) / comptonAngular(0.0, 511.0e3), 0.2, atol=0.02)
end