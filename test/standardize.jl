using Test
using NeXLCore

@testset "Standardize" begin
    @testset "Standardize KRatio" begin
        fe1 = KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
        ) 
        @test isstandard(fe1)
        fe2 = KRatio(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.3217, 0.02),
        )
        @test isstandard(fe2)
        ca1 = KRatio(
            characteristic(n"Ca", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Ca5(PO4)3F"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.7280, 0.02),
        )
        @test isstandard(ca1)
        ca2 = KRatio(
            characteristic(n"Ca", kbeta),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Ca5(PO4)3F"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.7326, 0.02),
        )
        @test isstandard(ca2)
        si1 = KRatio(
            characteristic(n"Si", ktransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"SiO2"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Si",
            uv(0.3801, 0.02),
        )
        @test isstandard(si1)
        @test NeXLCore.matches(KRatio(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
            ), fe2)
        @test !NeXLCore.matches(KRatio(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
            ), fe1)
        @test !NeXLCore.matches(KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.45, 0.02),
        ), ca1)
        stdized = standardize([ KRatio(
                characteristic(n"Fe", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"Fe",
                uv(0.45, 0.02)),
            KRatio(
                characteristic(n"Ca", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"CaF2",
                uv(0.32, 0.02)),
            KRatio(
                characteristic(n"F", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"LiF",
                uv(0.333, 0.02)) ],
            [ fe1, fe2, ca1, ca2, si1 ])
        @test value(stdized[1]) ==0.45/0.6517
        @test isequal(stdized[1].standard, mat"Fe2O3")
        @test value(stdized[2])==0.32/0.7280
        @test isequal(stdized[2].standard, mat"Ca5(PO4)3F")
        @test isequal(stdized[3].standard, mat"LiF")
        @test value(stdized[3])==0.333
    end
    @testset "Standardize KRatios" begin
        fe1 = KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
        ) 
        fe2 = KRatio(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.3217, 0.02),
        )
        ca1 = KRatio(
            characteristic(n"Ca", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Ca5(PO4)3F"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.7280, 0.02),
        )
        ca2 = KRatio(
            characteristic(n"Ca", kbeta),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Ca5(PO4)3F"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"CaF2",
            uv(0.7326, 0.02),
        )
        si1 = KRatio(
            characteristic(n"Si", ktransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"SiO2"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Si",
            uv(0.3801, 0.02),
        )
        @test NeXLCore.matches(KRatios(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            [ uv(0.6517, 0.02), uv(0.6517, 0.02) ]
            ), fe2)
        @test !NeXLCore.matches(KRatios(
            characteristic(n"Fe", ltransitions),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            [ uv(0.6517, 0.02), uv(0.6517, 0.02) ]
            ), fe1)
        @test !NeXLCore.matches(KRatios(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            [ uv(0.45, 0.02), uv(0.45, 0.02) ]
        ), ca1)
        stdized = standardize([ KRatios(
                characteristic(n"Fe", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"Fe",
                [ uv(0.45, 0.02), uv(0.49, 0.03)]),
            KRatios(
                characteristic(n"Ca", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"CaF2",
                [ uv(0.32, 0.02), uv(0.32, 0.02)]),
            KRatios(
                characteristic(n"F", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"LiF",
                [ uv(0.333, 0.02), uv(0.333, 0.02) ]) ],
            [ fe1, fe2, ca1, ca2, si1 ])
        @test value(stdized[1].kratios[1]) == 0.45/0.6517
        @test isapprox(σ(stdized[1].kratios[1]),0.037,atol=0.001)
        @test value(stdized[1].kratios[2]) == 0.49/0.6517
        @test isapprox(σ(stdized[1].kratios[2]),0.051,atol=0.001)
        @test isequal(stdized[1].standard, mat"Fe2O3")
        @test value(stdized[2].kratios[1])==0.32/0.7280
        @test isequal(stdized[2].standard, mat"Ca5(PO4)3F")
        @test isequal(stdized[3].standard, mat"LiF")
        @test value(stdized[3].kratios[1])==0.333
    end
    @testset "Standardize KRatios - mixed" begin
        fe1 = KRatio(
            characteristic(n"Fe", kalpha),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0), :Composition => mat"Fe2O3"),
            Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
            mat"Fe",
            uv(0.6517, 0.02),
        ) 
        stdized = standardize([ 
            KRatios(
                characteristic(n"Fe", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"Fe",
                [ 0.45, 0.49 ]),  # Float64
            KRatios(
                characteristic(n"Fe", kalpha),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
                mat"Fe",
                [ uv(0.45, 0.02), uv(0.49, 0.03) ]) ],  # UncertainValue
            [ fe1 ])
        # Check mixed UncertainValue and Float64
        @test stdized[1].kratios[1] == 0.45/0.6517
        @test σ(stdized[1].kratios[1]) == 0.0
        @test stdized[1].kratios[2] == 0.49/0.6517
        @test value(stdized[2].kratios[1]) == 0.45/0.6517
        @test isapprox(σ(stdized[2].kratios[1]),0.037,atol=0.001)
        @test isequal(stdized[1].standard, mat"Fe2O3")
    end
end