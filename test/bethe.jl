using Test

@testset "Bethe" begin
    @test isapprox(dEds(Bethe, 15.0e3, n"Fe", 8.0),-8.018e7, atol=0.01e7)
    @test isapprox(dEds(JoyLuo, 15.0e3, n"Fe", 8.0), -8.048e7, atol=0.01e7)

    @test isapprox(dEds(Bethe, 5.0e3, n"Al", 8.0),-2.154e8, atol=0.01e7)
    @test isapprox(dEds(JoyLuo, 5.0e3, n"Al", 8.0), -2.1696e8, atol=0.01e7)

    @test isapprox(dEds(Bethe, 1.0e3, n"U", 18.0),-1.475e8, atol=0.01e7)
    @test isapprox(dEds(JoyLuo, 1.0e3, n"U", 18.0), -4.5961e8, atol=0.01e7)

    @test dEds(Bethe, 1.0e2, n"U", 18.0) > 0.0
    @test isapprox(dEds(JoyLuo, 1.0e2, n"U", 18.0), -7.1987e8, atol=0.01e7)

    @test isapprox(range(JoyLuo, parse(Material, "BaTiSi3O9", density=3.65), 20.0e3)*1.0e4, 3.550, atol=0.01)
    @test isapprox(range(Bethe, parse(Material, "BaTiSi3O9", density=3.65), 20.0e3, true, emin=200.0)*1.0e4, 3.5760, atol=0.01)
    @test isapprox(range(Kanaya1972, parse(Material, "BaTiSi3O9", density=3.65), 20.0e3)*1.0e4, 3.32, atol=0.01)
end
