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

    v1 = [ -0.2605, 22.8789, 25.695 ]
    v2 = [ -13.7883, 21.1245, 25.662 ]
    v3 = [-4.3441, -0.1964, 25.820 ]
    (ϕ, θ) = compute_tilt(v1, v2, v3)
    @test isapprox(ϕ, 0.38933, atol=0.00001)
    @test isapprox(θ, 118.2447, atol=0.0001)
    
    v1 = [ 9.0, 0.0, 0.0 ]
    v2 = [ -8.0, 0.0, 0.0 ]
    v3 = [ 0.0, -10.0*cos(deg2rad(1.0)), 10.0*sin(deg2rad(1.0)) ]
    (ϕ, θ) = compute_tilt(v1, v2, v3)
    @test isapprox(ϕ, 1.0, atol=0.00001)
    @test isapprox(θ, 90.0, atol=0.0001)
end