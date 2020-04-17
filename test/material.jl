using Test
using NeXLCore

@testset "Materials" begin
    @testset "K411" begin
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

        @test k411[n"Ca"] == 0.110563
        @test k411[n"Ce"] == 0.0
        @test massfraction(k411)[n"Ca"] == k411[n"Ca"]
        @test k411[n"Ca"] == 0.110563

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
    @testset "Parse" begin
        albite = parse(Material, "AlNaSi3O8")
        @test isapprox(analyticaltotal(albite), 1.0, atol = 1.0e-8)
        @test isapprox(albite[n"O"], 0.488116, atol = 0.00001)
        @test isapprox(albite[n"Al"], 0.102895, atol = 0.00001)
        @test isapprox(albite[n"Si"], 0.321316, atol = 0.00001)
        @test isapprox(albite[n"Na"], 0.0876726, atol = 0.00001)

        apatite = parse(Material, "Ca5(PO4)3F")
        @test isapprox(apatite[n"O"], 0.38071, atol = 0.00001)
        @test isapprox(apatite[n"F"], 0.0376719, atol = 0.00001)
        @test isapprox(apatite[n"P"], 0.184257, atol = 0.00001)
        @test isapprox(apatite[n"Ca"], 0.397361, atol = 0.00001)

        apatite = parse(Material, "Ca5(PO4)3(OH)")
        @test isapprox(apatite[n"O"], 0.41407, atol = 0.00001)
        @test isapprox(apatite[n"H"], 0.0020066, atol = 0.00001)
        @test isapprox(apatite[n"P"], 0.184987, atol = 0.00001)
        @test isapprox(apatite[n"Ca"], 0.398936, atol = 0.00001)

        apatite = parse(Material, "Ca5(PO4)3OH")
        @test isapprox(apatite[n"O"], 0.41407, atol = 0.00001)
        @test isapprox(apatite[n"H"], 0.0020066, atol = 0.00001)
        @test isapprox(apatite[n"P"], 0.184987, atol = 0.00001)
        @test isapprox(apatite[n"Ca"], 0.398936, atol = 0.00001)

        apatite = parse(Material, "Ca5(PO4)3OH")
        @test isapprox(apatite[n"O"], 0.41407, atol = 0.00001)
        @test isapprox(apatite[n"H"], 0.0020066, atol = 0.00001)
        @test isapprox(apatite[n"P"], 0.184987, atol = 0.00001)
        @test isapprox(apatite[n"Ca"], 0.398936, atol = 0.00001)

        quartz1 = parse(Material, "Si23O46")
        quartz2 = parse(Material, "SiO2")
        @test all(e -> isapprox(quartz1[e], quartz2[e], atol = 0.00001), keys(quartz1))
        @test isapprox(quartz1, quartz2, atol = 1.0e-6)

        ss = parse(Material, "0.8*Fe+0.15*Ni+0.04*Cr")
        @test ss[n"Fe"] == 0.80
        @test ss[n"Ni"] == 0.15
        @test ss[n"Cr"] == 0.04

        m = parse(Material, "0.6*SiO2+0.4*Al2O3")
        @test isapprox(m[n"Si"], 0.6 * mat"SiO2"[n"Si"], atol = 1.0e-6)
        @test isapprox(m[n"Al"], 0.4 * mat"Al2O3"[n"Al"], atol = 1.0e-6)
        @test isapprox(m[n"O"], 0.6 * mat"SiO2"[n"O"] + 0.4 * mat"Al2O3"[n"O"], atol = 1.0e-6)

        @test isapprox(parse(Material, "Ba(Al2Si2O8)", name = "Paracelsian"), mat"Ba(Al2Si2O8)", atol = 1.0e-8)
        @test isapprox(mat"Ba(Al2Si2O8)"[n"Ba"], 0.36576, atol = 0.00001)

        para = Material(
            "Paracelsian",
            Dict(z(n"O") => 0.340906, z(n"Ba") => 0.36576, z(n"Si") => 0.149607, z(n"Al") => 0.143727),
            3.5,
        )
        @test isapprox(para, mat"Ba(Al2Si2O8)",atol=1.0e-5)
        para = material(
            "Paracelsian",
            Dict(n"O" => 0.340906, n"Ba" => 0.36576, n"Si" => 0.149607, n"Al" => 0.143727),
            3.5,
        )
        @test isapprox(para, mat"Ba(Al2Si2O8)",atol=1.0e-5)

        para = material(
            "Paracelsian",
            n"O" => 0.340906, n"Ba" => 0.36576, n"Si" => 0.149607, n"Al" => 0.143727,
            density = 3.5, description = "Barium-rich feldspar"
        )
        @test isapprox(para, mat"Ba(Al2Si2O8)",atol=1.0e-5)
        # Define a couple of additional properties...
        para[:Color]="Colorless to white, light yellow tint"
        para[:System]=:Monoclinic
        @test para[:Description]=="Barium-rich feldspar"
        @test para[:Density]==3.5
        @test density(para)==3.5
        @test para[:Color]=="Colorless to white, light yellow tint"
        @test para[:System]==:Monoclinic
    end
end
