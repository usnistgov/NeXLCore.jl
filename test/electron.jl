using Test
using NeXLCore


@testset "Electron" begin

    @test isapprox(λₑ(1.0), 1.227e-7, atol=0.001e-7)
    @test isapprox(λₑ(100.0), 1.227e-8, atol=0.001e-8)
    @test isapprox(λₑ(1.0e4), 1.227e-9, atol=0.001e-9)
    @test isapprox(λₑ(5.0e3), 0.0173*1.0e-7, atol=0.0001*1.0e-7)
    @test isapprox(λₑ(10.0e3), 0.0122*1.0e-7, atol=0.0001*1.0e-7)
    @test isapprox(λₑ(200.0e3), 0.0027*1.0e-7, atol=0.0001*1.0e-7)

    @test isapprox(kₑ(1.0), 0.01*sqrt(2.0*9.109e-31*1.602e-19*1.0)/1.05457e-34,rtol=1.0e-3)
    @test isapprox(kₑ(10.0), 0.01*sqrt(2.0*9.109e-31*1.602e-19*10.0)/1.05457e-34,rtol=1.0e-3)
    @test isapprox(kₑ(20.0), 0.01*sqrt(2.0*9.109e-31*1.602e-19*20.0)/1.05457e-34,rtol=1.0e-3)
    @test isapprox(kₑ(200.0), 0.01*sqrt(2.0*9.109e-31*1.602e-19*200.0)/1.05457e-34,rtol=1.0e-3)
    @test isapprox(kₑ(2000.0), 0.01*sqrt(2.0*9.109e-31*1.602e-19*2000.0)/1.05457e-34,rtol=1.0e-3)

    @test isapprox(Rₐ(n"Po"), 0.1208*1.0e-8, atol=0.0001*1.0e-8)
    @test isapprox(Rₐ(n"H"), 0.529*1.0e-8, atol=0.0001*1.0e-8)
    @test isapprox(Rₐ(n"Fe"), 0.1786*1.0e-8, atol=0.0001*1.0e-8)

    @test isapprox(mₑ, 0.511e6, atol= 0.001e6)  # Must be in eV!!!!
    @test isapprox(qₑ, 1.602e-19, atol= 0.001e-19)  # C 
end