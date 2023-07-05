using Test


@testset "Backscatter - Elemental" begin
    @test isapprox(η(LoveScott1978η, n"C", 10.0e3), 0.08432227, atol = 1.0e-5)
    @test isapprox(η(LoveScott1978η, n"Fe", 30.0e3), 0.2804824, atol = 1.0e-5)
    @test isapprox(η(LoveScott1978η, n"Ce", 10.0e3), 0.44044972, atol = 1.0e-5)
    @test isapprox(η(LoveScott1978η, n"U", 12.0e3), 0.5104889, atol = 1.0e-5)

    @test η(Tomlin1963, n"H", 10.0e3) > 0.0
    @test isapprox(η(Tomlin1963, n"C", 10.0e3), 0.04862657, atol = 1.0e-5)
    @test isapprox(η(Tomlin1963, n"Fe", 30.0e3), 0.2930160, atol = 1.0e-5)
    @test isapprox(η(Tomlin1963, n"Ce", 10.0e3), 0.42674050, atol = 1.0e-5)
    @test isapprox(η(Tomlin1963, n"U", 12.0e3), 0.5036314, atol = 1.0e-5)

    @test all(map(z -> η(August1989η, elements[z], 15.0e3) == η(elements[z], 15.0e3), 1:99))
end

@testset "Backscatter - Material" begin

    @test isapprox(
        η(LoveScott1978η, n"C", 10.0e3),
        η(LoveScott1978η, mat"0.9999*C+0.0001*H", 10.0e3),
        atol = 1.0e-5,
    )
    @test isapprox(
        η(LoveScott1978η, n"U", 30.0e3),
        η(LoveScott1978η, mat"0.99999*U+0.00001*H", 30.0e3),
        atol = 1.0e-5,
    )
    @test isapprox(
        η(LoveScott1978η, n"C", 10.0e3),
        η(LoveScott1978η, mat"0.9999*C+0.0001*H", 10.0e3),
        atol = 1.0e-5,
    )

    @test isapprox(
        NeXLCore.electronfraction(n"Au", mat"0.756*Au+0.244*Cu"),
        0.731,
        atol = 0.001,
    )
    @test isapprox(NeXLCore.electronfraction(n"Th", mat"ThSiO4"), 0.6618, atol = 0.0001)
    @test isapprox(NeXLCore.electronfraction(n"Si", mat"ThSiO4"), 0.1029, atol = 0.0001)

    @test isapprox(η(LoveScott1978η, mat"0.224*Au+0.775*Ag", 20.0e3), 0.4271, atol = 1.0e-3)
    @test isapprox(η(LoveScott1978η, mat"0.80*Au+0.199*Ag", 20.0e3), 0.4761, atol = 1.0e-3)

    @test isapprox(η(LoveScott1978η, mat"0.201*Au+0.798*Cu", 20.0e3), 0.3511, atol = 1.0e-3)
    @test isapprox(η(LoveScott1978η, mat"0.801*Au+0.198*Cu", 20.0e3), 0.4610, atol = 1.0e-3)

    @test isapprox(zbar(mat"PbS), 66.053963, atol=0.0001)
    @test isapprox(zbar(mat"ZrSiO4), 20.5635, atol=0.0001)
end
