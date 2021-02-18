using NeXLCore
using Test

@testset "Staging" begin
    stg_coords = Dict(:X=>10.0, :Y=>5.0)
    img_coords = Dict(:X=>0.027777, :Y=>-0.016666 )
    rot = deg2rad(50.0)
    i2s = image2stage(NeXLCore.DefaultStageMapping, stg_coords, img_coords, rot)
    s2i = stage2image(NeXLCore.DefaultStageMapping, stg_coords, i2s, rot)
    @test isapprox(s2i[:X], img_coords[:X], atol=1.0e-6)    
    @test isapprox(s2i[:Y], img_coords[:Y], atol=1.0e-6)
end