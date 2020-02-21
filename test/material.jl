using Test
using NeXLCore

@testset "Materials" begin

        k411 = NeXLCore.Material(
                "K411",
                Dict(
                        z(n"O") => 0.423667,
                        z(n"Mg") => 0.0884654,
                        z(n"Si") => 0.253817,
                        z(n"Ca") => 0.110563,
                        z(n"Fe") => 0.112087,
                ),
                3.5,
        )

        @test massfraction(n"Ca", k411) == 0.110563
        @test massfraction(n"Ce", k411) == 0.0
        @test massfraction(k411)[n"Ca"] == massfraction(n"Ca", k411)
        @test massfraction(n"Ca", k411) == 0.110563

        nk411 = normalizedmassfraction(k411)

        @test isapprox(nk411[n"O"], 0.428553, atol = 1.0e-6)
        @test isapprox(nk411[n"Mg"], 0.0894855, atol = 1.0e-6)
        @test isapprox(nk411[n"Si"], 0.256744, atol = 1.0e-6)
        @test isapprox(nk411[n"Ca"], 0.111838, atol = 1.0e-6)
        @test isapprox(nk411[n"Fe"], 0.11338, atol = 1.0e-6)

        @test analyticaltotal(k411) == 0.9885994

        ak411 = atomicfraction(k411)
        @test isapprox(ak411[n"O"], 0.602876, atol = 1.0e-5)
        @test isapprox(ak411[n"Mg"], 0.0828675, atol = 1.0e-5)
        @test isapprox(ak411[n"Si"], 0.205753, atol = 1.0e-5)
        @test isapprox(ak411[n"Ca"], 0.0628072, atol = 1.0e-5)
        @test isapprox(ak411[n"Fe"], 0.0456961, atol = 1.0e-5)

        k411_2 = material(
                "K411",
                Dict(n"O" => 0.423667, n"Mg" => 0.0884654, n"Si" => 0.253817, n"Ca" => 0.110563, n"Fe" => 0.112087),
                3.5,
                Dict(n"Ca" => 41.0, n"O" => 16.1),
        )

        @test a(n"Ca", k411) == a(n"Ca")
        @test a(n"Ca", k411_2) == 41.0
        @test a(n"O", k411_2) == 16.1
        @test a(n"Si", k411) == a(n"Si")
        @test a(n"Fe", k411) == a(n"Fe")
        @test a(n"Si", k411_2) == a(n"Si")
        @test a(n"Fe", k411_2) == a(n"Fe")

end
