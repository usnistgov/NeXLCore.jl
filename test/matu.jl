using Test
using NeXLCore

@testset "MatUncertain" begin
    mat = "Au60Ag40"
    inplbls = [
        MassFractionLabel(mat, n"Ag"),
        AtomicWeightLabel(mat, n"Ag"),
        MassFractionLabel(mat, n"Au"),
        AtomicWeightLabel(mat, n"Au"),
    ]
    inpvals = Float64[0.4020, 107.87, 0.5950, 196.97]
    inpcovs = Float64[(0.0090)^2 0 0 0; 0 0 0 0; 0 0 (0.0120)^2 0; 0 0 0 0]
    inputs = uvs(inplbls, inpvals, inpcovs)

    @testset "AuAg - MFtoAF" begin

        afs = propagate(MFtoAF(mat), inputs)

        @test isapprox(value(afs[AtomicFractionLabel(mat, n"Ag")]), 0.5523, atol = 0.0001)
        @test isapprox(value(afs[AtomicFractionLabel(mat, n"Au")]), 0.4477, atol = 0.0001)

        @test isapprox(σ(afs[AtomicFractionLabel(mat, n"Ag")]), 0.0075, atol = 0.0001)
        @test isapprox(σ(afs[AtomicFractionLabel(mat, n"Au")]), 0.0075, atol = 0.0001)
        @test isapprox(
            NeXLUncertainties.correlation(
                afs,
                AtomicFractionLabel(mat, n"Ag"),
                AtomicFractionLabel(mat, n"Au"),
                
            ),
            -1.0,
            atol = 0.0001,
        )
    end

    @testset "AuAg - MFtoNMF" begin

        afs = propagate(MFtoNMF(mat), inputs)

        @test isapprox(value(afs[NormMassFractionLabel(mat, n"Ag")]), 0.4032, atol = 0.0001)
        @test isapprox(value(afs[NormMassFractionLabel(mat, n"Au")]), 0.5968, atol = 0.0001)

        @test isapprox(σ(afs[NormMassFractionLabel(mat, n"Ag")]), 0.0073, atol = 0.0001)
        @test isapprox(σ(afs[NormMassFractionLabel(mat, n"Au")]), 0.0073, atol = 0.0001)
        @test isapprox(
            NeXLUncertainties.correlation(
                afs,
                NormMassFractionLabel(mat, n"Ag"),
                NormMassFractionLabel(mat, n"Au"),
            ),
            -1.0,
            atol = 0.0001,
        )
    end

    @testset "AuAg - MatStats" begin

        afs = propagate(MatStats(mat), inputs)

        @test isapprox(value(afs[MeanZ(mat)]), 65.8990, atol = 0.0001)
        @test isapprox(value(afs[MeanAz(mat)]), 160.56, atol = 0.01)

        @test isapprox(σ(afs[MeanZ(mat)]), 1.0381, atol = 0.0001)
        @test isapprox(σ(afs[MeanAz(mat)]), 2.5552, atol = 0.0001)
        @test isapprox(
            NeXLUncertainties.correlation(afs, MeanZ(mat), MeanAz(mat)),
            1.0,
            atol = 0.001,
        )
    end

    @testset "AuAg - Comp" begin
        afs = mf2comp(mat, inputs)

        @test isapprox(value(afs[AtomicFractionLabel(mat, n"Ag")]), 0.5523, atol = 0.0001)
        @test isapprox(value(afs[AtomicFractionLabel(mat, n"Au")]), 0.4477, atol = 0.0001)

        @test isapprox(σ(afs[AtomicFractionLabel(mat, n"Ag")]), 0.0075, atol = 0.0001)
        @test isapprox(σ(afs[AtomicFractionLabel(mat, n"Au")]), 0.0075, atol = 0.0001)
        @test isapprox(
            NeXLUncertainties.correlation(
                afs,
                AtomicFractionLabel(mat, n"Ag"),
                AtomicFractionLabel(mat, n"Au")
            ),
            -1.0,
            atol = 0.0001,
        )

        @test isapprox(value(afs[NormMassFractionLabel(mat, n"Ag")]), 0.4032, atol = 0.0001)
        @test isapprox(value(afs[NormMassFractionLabel(mat, n"Au")]), 0.5968, atol = 0.0001)

        @test isapprox(σ(afs[NormMassFractionLabel(mat, n"Ag")]), 0.0073, atol = 0.0001)
        @test isapprox(σ(afs[NormMassFractionLabel(mat, n"Au")]), 0.0073, atol = 0.0001)
        @test isapprox(
            NeXLUncertainties.correlation(
                afs,
                NormMassFractionLabel(mat, n"Ag"),
                NormMassFractionLabel(mat, n"Au")
            ),
            -1.0,
            atol = 0.0001,
        )

        @test isapprox(value(afs[MeanZ(mat)]), 65.8990, atol = 0.0001)
        @test isapprox(value(afs[MeanAz(mat)]), 160.56, atol = 0.01)

        @test isapprox(σ(afs[MeanZ(mat)]), 1.0381, atol = 0.0001)
        @test isapprox(σ(afs[MeanAz(mat)]), 2.5552, atol = 0.0001)
        @test isapprox(
            NeXLUncertainties.correlation(afs, MeanZ(mat), MeanAz(mat)),
            1.0,
            atol = 0.001,
        )
    end

    @testset "AuAg - MacU" begin

        mat, cxr = material("Au60Ag40", n"Ag" => 0.4020, n"Au" => 0.5950), n"O K-L3"
        inplbls = [
            MassFractionLabel(mat.name, n"Ag"),
            AtomicWeightLabel(mat.name, n"Ag"),
            MassFractionLabel(mat.name, n"Au"),
            AtomicWeightLabel(mat.name, n"Au"),
            μoρElementLabel(n"Ag", cxr),
            μoρElementLabel(n"Au", cxr),
        ]
        macAg, macAu = macU(n"Ag", cxr), macU(n"Au", cxr)
        @test isapprox(value(macAg), mac(n"Ag", cxr), rtol = 1.0e-6)
        @test isapprox(value(macAu), mac(n"Au", cxr), rtol = 1.0e-6)

        inpvals = Float64[0.4020, 107.87, 0.5950, 196.97, value(macAg), value(macAu)]
        inpcovs = Float64[
            (0.0090)^2 0 0 0 0 0
            0 0 0 0 0 0
            0 0 (0.0120)^2 0 0 0
            0 0 0 0 0 0
            0 0 0 0 variance(macAg) 0
            0 0 0 0 0 variance(macAu)
        ]

        inputs = uvs(inplbls, inpvals, inpcovs)
        model = μoρMaterial(mat.name, cxr)
        res = model(inputs)
        mc_res = mcpropagate(
            model,
            inputs,
            1000,
            parallel = false,
            rng = MersenneTwister(0xBADF00D),
        )
        #println(res)
        #println(mc_res)
        lbl = μoρLabel(mat.name, cxr)
        @test isapprox(value(res[lbl]), mac(mat, cxr), rtol = 1.0e-6)
        @test isapprox(value(res[lbl]), value(mc_res[lbl]), rtol = 0.01)
        @test isapprox(σ(res[lbl]), σ(mc_res[lbl]), rtol = 0.01)
    end
end
