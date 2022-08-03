using Test
using NeXLCore
using Unitful

@testset "Electron" begin

    @test isapprox(λₑ(1.0), 1.227e-7, atol = 0.001e-7)
    @test isapprox(λₑ(100.0), 1.227e-8, atol = 0.001e-8)
    @test isapprox(λₑ(1.0e4), 1.22e-9, atol = 0.001e-9)
    @test isapprox(λₑ(5.0e3), 0.0173 * 1.0e-7, atol = 0.0001 * 1.0e-7)
    @test isapprox(λₑ(10.0e3), 0.01227 * 1.0e-7, atol = 0.0001 * 1.0e-7)
    @test isapprox(λₑ(200.0e3), 0.0025 * 1.0e-7, atol = 0.0001 * 1.0e-7)
    @test isapprox(λₑ(51.0u"GeV"), 2.4e-17u"m", atol=0.1e-17u"m")

    @test isapprox(
        kₑ(1.0),
        0.01 * sqrt(2.0 * 9.109e-31 * 1.602e-19 * 1.0) / 1.05457e-34,
        rtol = 1.0e-3,
    )
    @test isapprox(
        kₑ(10.0),
        0.01 * sqrt(2.0 * 9.109e-31 * 1.602e-19 * 10.0) / 1.05457e-34,
        rtol = 1.0e-3,
    )
    @test isapprox(
        kₑ(20.0),
        0.01 * sqrt(2.0 * 9.109e-31 * 1.602e-19 * 20.0) / 1.05457e-34,
        rtol = 1.0e-3,
    )
    @test isapprox(
        kₑ(200.0),
        0.01 * sqrt(2.0 * 9.109e-31 * 1.602e-19 * 200.0) / 1.05457e-34,
        rtol = 1.0e-3,
    )
    @test isapprox(kₑ(2000.0), 2.2933e9, rtol = 1.0e-3)
    @test isapprox(mₑ, 0.511e6, atol = 0.001e6)  # Must be in eV!!!!

    @test isapprox(NeXLCore.electrons_per_second(1.0u"nA"), 6.24e9*u"s^-1", atol=0.01e9*u"s^-1")
    @test isapprox(NeXLCore.electrons_per_second(1.0), 6.24e9, atol=0.01e9)

    @test isapprox(NeXLCore.vₑ(0.511e6), 2.596e10, atol = 0.001e10)
    @test isapprox(NeXLCore.vₑ(0.511e6u"eV"), 2.596e10u"cm/s", atol = 0.001e10u"cm/s")

    @test isapprox(NeXLCore.vₑ(51e9), ustrip(SpeedOfLightInVacuum |> u"cm/s"), rtol = 0.001)
    @test isapprox(NeXLCore.vₑ(51u"GeV"), SpeedOfLightInVacuum, rtol = 0.001)

    @test isapprox(NeXLCore.γₑ(0.999*SpeedOfLightInVacuum), 22.36, rtol = 0.001)
    @test isapprox(NeXLCore.γₑ(0.01*SpeedOfLightInVacuum), 1.0, rtol = 0.001)
end