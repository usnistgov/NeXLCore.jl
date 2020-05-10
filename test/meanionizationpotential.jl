using Test
using NeXLCore

@testset "MIP" begin
    @test all(map(z->isapprox(J(Bloch1933, z)/z,13.5,rtol=1.0e-8),1:99))
    @test all(map(z->isapprox(J(Jensen1937, z)/z,9.0*(1.0+0.5*z^(-2.0/3.0)),rtol=1.0e-8),1:99))
    @test all(map(z->isapprox(J(Wilson1941, z)/z,11.5,rtol=1.0e-8),1:99))
    @test all(map(z->isapprox(J(Sternheimer1964, z)/z,9.76+58.82*z^-1.19,rtol=1.0e-8),12:99))
    @test all(map(z->isapprox(J(Springer1967, z)/z,9.0*(1.0+z^(-2.0/3.0))+0.03*z,rtol=1.0e-8),1:99))
    @test all(map(z->isapprox(J(Zeller1975, z)/z, 10.04 + 8.25*exp(-z/11.22),rtol=1.0e-8),1:99))
    @test all(map(z->isapprox(J(Brizuela1990, z)/z,22.4*z^-0.172,rtol=1.0e-8),1:99))
    @test J(Berger1983, n"C")==78.0
    @test J(Berger1983, n"Fe")==286.0
    @test J(Berger1983, n"Pb")==823.0
    @test J(Berger1983, n"Pu")==921.0
end
