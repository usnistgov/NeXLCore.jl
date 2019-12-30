using Test
using NeXLCore

@testset "MatUncertain" begin
    mat = "Au60Ag40"
    inplbls = [ MassFractionLabel(mat,n"Ag"), AtomicWeightLabel(mat,n"Ag"),
                MassFractionLabel(mat,n"Au"), AtomicWeightLabel(mat,n"Au") ]
    inpvals = Float64[ 0.4020, 107.87, 0.5950, 196.97 ]
    inpcovs = Float64[ (0.0090)^2 0 0 0; 0 0 0 0; 0 0 (0.0120)^2 0; 0 0 0 0 ]
    inputs = uvs(inplbls, inpvals, inpcovs)

    @testset "AuAg - MFtoAF" begin

        afs = propagate(MFtoAF(mat), inputs)

        @test isapprox(value(afs[AtomicFractionLabel(mat,n"Ag")]),0.5523,atol=0.0001)
        @test isapprox(value(afs[AtomicFractionLabel(mat,n"Au")]),0.4477,atol=0.0001)

        @test isapprox(σ(afs[AtomicFractionLabel(mat,n"Ag")]),0.0075,atol=0.0001)
        @test isapprox(σ(afs[AtomicFractionLabel(mat,n"Au")]),0.0075,atol=0.0001)
        @test isapprox(NeXLUncertainties.correlation(AtomicFractionLabel(mat,n"Ag"), AtomicFractionLabel(mat,n"Au"), afs), -1.0,atol=0.0001)
    end;

    @testset "AuAg - MFtoNMF" begin

        afs = propagate(MFtoNMF(mat), inputs)

        @test isapprox(value(afs[NormMassFractionLabel(mat,n"Ag")]),0.4032,atol=0.0001)
        @test isapprox(value(afs[NormMassFractionLabel(mat,n"Au")]),0.5968,atol=0.0001)

        @test isapprox(σ(afs[NormMassFractionLabel(mat,n"Ag")]),0.0073,atol=0.0001)
        @test isapprox(σ(afs[NormMassFractionLabel(mat,n"Au")]),0.0073,atol=0.0001)
        @test isapprox(NeXLUncertainties.correlation(NormMassFractionLabel(mat,n"Ag"), NormMassFractionLabel(mat,n"Au"), afs), -1.0,atol=0.0001)
    end;

    @testset "AuAg - MatStats" begin

        afs = propagate(MatStats(mat), inputs)

        @test isapprox(value(afs[MeanZ(mat)]), 65.8990, atol=0.0001)
        @test isapprox(value(afs[MeanAz(mat)]), 160.56, atol=0.01)

        @test isapprox(σ(afs[MeanZ(mat)]), 1.0381, atol=0.0001)
        @test isapprox(σ(afs[MeanAz(mat)]), 2.5552, atol=0.0001)
        @test isapprox(NeXLUncertainties.correlation(MeanZ(mat), MeanAz(mat), afs), 1.0, atol=0.001)
    end;

    @testset "AuAg - Comp" begin
        afs = mf2comp( mat, inputs)

        @test isapprox(value(afs[AtomicFractionLabel(mat,n"Ag")]),0.5523,atol=0.0001)
        @test isapprox(value(afs[AtomicFractionLabel(mat,n"Au")]),0.4477,atol=0.0001)

        @test isapprox(σ(afs[AtomicFractionLabel(mat,n"Ag")]),0.0075,atol=0.0001)
        @test isapprox(σ(afs[AtomicFractionLabel(mat,n"Au")]),0.0075,atol=0.0001)
        @test isapprox(NeXLUncertainties.correlation(AtomicFractionLabel(mat,n"Ag"), AtomicFractionLabel(mat,n"Au"), afs), -1.0,atol=0.0001)

        @test isapprox(value(afs[NormMassFractionLabel(mat,n"Ag")]),0.4032,atol=0.0001)
        @test isapprox(value(afs[NormMassFractionLabel(mat,n"Au")]),0.5968,atol=0.0001)

        @test isapprox(σ(afs[NormMassFractionLabel(mat,n"Ag")]),0.0073,atol=0.0001)
        @test isapprox(σ(afs[NormMassFractionLabel(mat,n"Au")]),0.0073,atol=0.0001)
        @test isapprox(NeXLUncertainties.correlation(NormMassFractionLabel(mat,n"Ag"), NormMassFractionLabel(mat,n"Au"), afs), -1.0,atol=0.0001)

        @test isapprox(value(afs[MeanZ(mat)]), 65.8990, atol=0.0001)
        @test isapprox(value(afs[MeanAz(mat)]), 160.56, atol=0.01)

        @test isapprox(σ(afs[MeanZ(mat)]), 1.0381, atol=0.0001)
        @test isapprox(σ(afs[MeanAz(mat)]), 2.5552, atol=0.0001)
        @test isapprox(NeXLUncertainties.correlation(MeanZ(mat), MeanAz(mat), afs), 1.0, atol=0.001)
    end;

end;
