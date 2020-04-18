using Test
using NeXLCore

@testset "KRatio" begin
    kr = KRatio( [ n"Fe K-L2", n"Fe K-L3" ],
            Dict(:BeamEnergy=>20.0e3,:TakeOffAngle=>deg2rad(40.0)),
            Dict(:BeamEnergy=>15.0e3,:TakeOffAngle=>deg2rad(38.0)),
            mat"Fe2O3",
            uv(0.77,0.01))
    @test n"Fe K-L2" in kr.lines
    @test !(n"Ca K-L2" in kr.lines)
    @test kr.unkProps[:BeamEnergy]==20.0e3
    @test kr.unkProps[:TakeOffAngle]==deg2rad(40.0)
    @test kr.stdProps[:BeamEnergy]==15.0e3
    @test kr.stdProps[:TakeOffAngle]==deg2rad(38.0)
    @test isapprox(kr.standard[n"Fe"],0.6994315143211048,atol=1e-8)
    @test value(kr.kratio)==0.77
    @test σ(kr.kratio)==0.01
    @test value(nonnegk(kr))==0.77
    @test σ(nonnegk(kr))==0.01
    kr = KRatio( [ n"Fe K-L2", n"Fe K-L3" ],
            Dict(:BeamEnergy=>20.0e3,:TakeOffAngle=>deg2rad(40.0)),
            Dict(:BeamEnergy=>15.0e3,:TakeOffAngle=>deg2rad(38.0)),
            mat"Fe2O3",
            uv(-0.01,0.01))
    @test value(nonnegk(kr))==0.0
    @test σ(nonnegk(kr))==0.01
    @test element(kr)==n"Fe"
    kr2 = KRatio( [ n"O K-L3" ],
            Dict(:BeamEnergy=>20.0e3,:TakeOffAngle=>deg2rad(40.0)),
            Dict(:BeamEnergy=>15.0e3,:TakeOffAngle=>deg2rad(38.0)),
            mat"Fe2O3",
            uv(0.6,0.01))
    krs=[kr, kr2]
    @test length(elms(krs))==2
    @test n"Fe" in elms(krs)
    @test n"O" in elms(krs)
    @test !(n"C" in elms(krs))
    @test_throws ErrorException KRatio( [ n"O K-L3", n"F K-L3" ],
            Dict(:BeamEnergy=>20.0e3,:TakeOffAngle=>deg2rad(40.0)),
            Dict(:BeamEnergy=>15.0e3,:TakeOffAngle=>deg2rad(38.0)),
            mat"Fe2O3",
            uv(0.6,0.01))
    @test_throws ErrorException KRatio( [ n"F K-L3" ],
            Dict(:BeamEnergy=>20.0e3,:TakeOffAngle=>deg2rad(40.0)),
            Dict(:BeamEnergy=>15.0e3,:TakeOffAngle=>deg2rad(38.0)),
            mat"Fe2O3",
            uv(0.6,0.01))
end
