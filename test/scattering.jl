using Test
using NeXLCore

@testset "Liljequist1989 tests" begin
    # Almost Table II
    @test isapprox(σₜ(Liljequist1989, n"Al", 100.0), 14.001*NeXLCore.a₀^2, atol = 0.001) # 14.0
    @test isapprox(σₜ(Liljequist1989, n"Al", 500.0), 5.879*NeXLCore.a₀^2, atol = 0.001) # 5.89
    @test isapprox(σₜ(Liljequist1989, n"Al", 1000.0), 3.523*NeXLCore.a₀^2, atol = 0.001) # 3.53
    @test isapprox(σₜ(Liljequist1989, n"Al", 5000.0), 0.894*NeXLCore.a₀^2, atol=0.001) # 0.899
    @test isapprox(σₜ(Liljequist1989, n"Al", 10000.0), 0.4695*NeXLCore.a₀^2, atol = 0.001) # 0.474
    @test isapprox(σₜ(Liljequist1989, n"Al", 50000.0), 0.1037*NeXLCore.a₀^2, atol = 0.001) # 0.109
    @test isapprox(σₜ(Liljequist1989, n"Al", 100000.0), 0.0560*NeXLCore.a₀^2, atol = 0.001) # 0.0615

    @test isapprox(σₜ(Liljequist1989, n"Au", 100.0), 6.775*NeXLCore.a₀^2, atol = 0.001) # 6.78
    @test isapprox(σₜ(Liljequist1989, n"Au", 500.0), 8.738*NeXLCore.a₀^2, atol = 0.001) # 8.75
    @test isapprox(σₜ(Liljequist1989, n"Au", 1000.0), 9.303*NeXLCore.a₀^2, atol = 0.001) # 9.32
    @test isapprox(σₜ(Liljequist1989, n"Au", 5000.0), 4.762*NeXLCore.a₀^2, atol=0.001) # 4.79
    @test isapprox(σₜ(Liljequist1989, n"Au", 10000.0), 3.095*NeXLCore.a₀^2, atol = 0.001) # 3.13
    @test isapprox(σₜ(Liljequist1989, n"Au", 50000.0), 1.041*NeXLCore.a₀^2, atol = 0.001) # 1.08
    @test isapprox(σₜ(Liljequist1989, n"Au", 100000.0), 0.630*NeXLCore.a₀^2, atol = 0.001) # 0.692
end