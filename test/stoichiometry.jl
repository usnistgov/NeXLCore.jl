using Test

using NeXLCore

@testset "Stoichiometry" begin

    # Check against DTSA-II
    @test isapprox(
        obystoichiometry(
            n"Mg" => 0.116567,
            n"Al" => 0.0490615,
            n"Si" => 0.211982,
            n"Ca" => 0.10899,
            n"Fe" => 0.0774196,
        ),
        0.4276,
        atol = 0.0001,
    )

    sio2 = asoxide(n"Si")
    @test name(sio2) == "SiO₂"
    asio2 = atomicfraction("SiO2", n"Si" => 1, n"O" => 2)
    asio2 = atomicfraction("SiO₂", n"Si" => 1, n"O" => 2)
    @test isapprox(sio2[n"Si"], asio2[n"Si"], atol = 1.0e-6)
    @test isapprox(sio2[n"O"], asio2[n"O"], atol = 1.0e-6)

    mix = sum(asoxide(           
        n"Mg" => 0.116567,
        n"Al" => 0.0490615,
        n"Si" => 0.211982,
        n"Ca" => 0.10899,
        n"Fe" => 0.0774196))
    
    @test isapprox(mix[n"O"], 0.4276, atol = 0.0001)
    @test isapprox(mix[n"Al"], 0.0490615, atol = 0.0001)
    @test isapprox(mix[n"Fe"], 0.0774196, atol = 0.0001)


    ox = sum(asoxide(filter(elm->elm[1]!=n"O", massfraction(srm470_k412))), name="K412")
    @test all(isapprox(ox[elm], srm470_k412[elm], atol=1.0e-8) for elm in keys(srm470_k412))

    ox = sum(asoxide(filter(elm->elm[1]!=n"O", massfraction(srm470_k411))), name="K411")
    @test all(isapprox(ox[elm], srm470_k411[elm], atol=1.0e-8) for elm in keys(srm470_k411))

    m=sum(asoxide(n"Al"=>0.1029, n"Na"=>0.0877, n"Si"=>0.3213), name="Albite", density=5.0, description="a common mineral", pedigree="Natural")
    @test isapprox(m[n"O"], 0.4881, atol=0.0001)
    @test isapprox(m[n"Si"], 0.3212, atol=0.0001)
    @test isapprox(m[n"Na"], 0.0877, atol=0.0001)
    @test m[:Density] == 5.0
    @test m[:Description] == "a common mineral"
    @test m[:Pedigree] == "Natural"
    @test name(m) == "Albite"

    amix =
        0.4 * atomicfraction("SiO2", n"Si" => 1, n"O" => 2) +
        0.6 * atomicfraction("Al2O3", n"Al" => 2, n"O" => 3)
    @test isapprox(amix[n"O"], 0.4955, atol = 0.0001)
    @test isapprox(amix[n"Al"], 0.3176, atol = 0.0001)
    @test isapprox(amix[n"Si"], 0.1870, atol = 0.0001)

    # Check against DTSA-II
    @test isapprox(amix[n"Si"], 0.1870, atol = 0.0001)
    @test isapprox(amix[n"Al"], 0.3176, atol = 0.0001)
    @test isapprox(amix[n"O"], 0.4955, atol = 0.0001)
end
