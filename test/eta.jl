using Test


@testset "Backscatter" begin
    @test isapprox(η(LoveScott1978η, n"C", 10.0e3),0.08432227,atol=1.0e-5)
    @test isapprox(η(LoveScott1978η, n"Fe", 30.0e3),0.2804824,atol=1.0e-5)
    @test isapprox(η(LoveScott1978η, n"Ce", 10.0e3),0.44044972,atol=1.0e-5)
    @test isapprox(η(LoveScott1978η, n"U", 12.0e3),0.5104889,atol=1.0e-5)

    @test isapprox(η(Tomlin1963, n"C", 10.0e3),0.08432227,atol=1.0e-5)
    @test isapprox(η(Tomlin1963, n"Fe", 30.0e3),0.2804824,atol=1.0e-5)
    @test isapprox(η(Tomlin1963, n"Ce", 10.0e3),0.44044972,atol=1.0e-5)
    @test isapprox(η(Tomlin1963, n"U", 12.0e3),0.5104889,atol=1.0e-5)


end
