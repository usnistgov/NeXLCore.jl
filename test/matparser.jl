using Test

@testset "Material Parser" begin
    @test NeXLCore._mp_level5("SiO2")  == Dict{Element, Int64}(n"Si" => 1, n"O" => 2)
    @test NeXLCore._mp_level5("Al2O3") == Dict{Element, Int64}(n"Al" => 2, n"O" => 3)
    @test NeXLCore._mp_level5("NaAlSi3O8") == Dict{Element, Int64}(n"Na" => 1, n"Al" => 1, n"Si" => 3, n"O" => 8)
    @test NeXLCore._mp_level5("") == Dict{Element,Int64}()
    @test_throws ErrorException NeXLCore._mp_level5("SiLO2")
    @test_throws ErrorException NeXLCore._mp_level5("SiO2_")
    @test_throws ErrorException NeXLCore._mp_level5("Si_O2")
    @test_throws ErrorException NeXLCore._mp_level5("9SiO2")

    @test NeXLCore._mp_level4("Ca5(PO4)3F") == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"F" => 1, n"O" => 12)
    @test NeXLCore._mp_level4("Ca5(PO4)3OH") == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"H" => 1, n"O" => 13)
    @test NeXLCore._mp_level4("(Ca)5(PO4)3OH") == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"H" => 1, n"O" => 13)
    @test NeXLCore._mp_level4("(Ca5)(PO4)3(OH)")  == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"H" => 1, n"O" => 13)
    @test NeXLCore._mp_level4("(Ca5)(PO4)3OH")  == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"H" => 1, n"O" => 13)
    @test_throws ErrorException NeXLCore._mp_level4("(Ca5(PO4)3OH")
    @test_throws ErrorException NeXLCore._mp_level4("Ca5((PO4)3OH")
    @test_throws ErrorException NeXLCore._mp_level4("Ca5)(PO4)3OH")

    @test NeXLCore._mp_level3("Ca5(PO4)3⋅OH") == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"H" => 1, n"O" => 13)
    @test NeXLCore._mp_level3("Ca5(PO4)3⋅1OH") == Dict{Element, Int}(n"P" => 3, n"Ca" => 5, n"H" => 1, n"O" => 13)

    function cda(d1, d2)
        all(isapprox(uv(get(d1,el, zero(valtype(d1)))), uv(get(d2,el,zero(valtype(d2)))), atol=1.0e-4) for el in intersect(keys(d1), keys(d2)))
    end

    # Check a couple of materials (Compare to DTSA-II)
    @test cda(NeXLCore.mmm_arsenopyrite.massfraction, Dict(n"S"=>0.1928, n"Fe"=>0.3416,n"Co"=>0.0002,n"As"=>0.4610))
    @test cda(NeXLCore.mmm_biotite.massfraction, Dict(n"H"=>0.0043,n"O"=>0.4372,n"Mg"=>0.1186,n"Al"=>0.0846,n"Si"=>0.1791,n"K"=>0.0835,n"Ti"=>0.0094,n"Fe"=>0.0780))
    @test cda(Dict(elm=>value(v) for (elm, v) in NeXLCore.srm470_k412.massfraction), Dict(n"O"=>0.42758,n"Mg"=>0.116567,n"Al"=>0.0490615,n"Si"=>0.211982,n"Ca"=>0.10899,n"Fe"=>0.077419))
    
    @test cda(NeXLCore._mp_level1("Ca5(PO4)3⋅1OH", Dict{Element,Float64}(), s->nothing), Dict{Element, Float64}(n"P" => 0.184989, n"Ca" => 0.398942, n"H" => 0.00200674, n"O" => 0.414062))

    # Test the lookup function for named materials
    function luf(s)
        d = Dict(
            "K412"=>NeXLCore.srm470_k412.massfraction, 
            "K411"=>NeXLCore.srm470_k411.massfraction,
            "Arsenopyrite"=>NeXLCore.mmm_arsenopyrite.massfraction
        ) 
        return get(d, s, nothing)
    end
    res1 = Dict(elm=>0.8*val for (elm,val) in srm470_k412.massfraction)
    res2 = Dict(elm=>0.2*val for (elm,val) in srm470_k411.massfraction)
    res = Dict(elm=>sum(uv(get(res1, elm, zero(valtype(res1)))), uv(get(res2, elm, zero(valtype(res2))))) for elm in union(keys(res1),keys(res2)))
    @test cda(parse(Material, "0.8*K412+0.2*K411", lookup=luf).massfraction, res)

    res1 = Dict(elm=>0.34*val for (elm,val) in NeXLCore.mmm_arsenopyrite.massfraction)
    res2 = Dict(elm=>0.62*val for (elm,val) in mat"NaAlSi3O8".massfraction)
    res = Dict(elm=>sum(uv(get(res1, elm, zero(valtype(res1)))), uv(get(res2, elm, zero(valtype(res2))))) for elm in union(keys(res1),keys(res2)))
    @test cda(parse(Material, "0.34*Arsenopyrite+0.62*NaAlSi3O8", lookup=luf).massfraction, res)

    @test isapprox(mat"RbTiO.PO4"[n"P"], 0.1268, atol=0.0001)
    @test isapprox(mat"Ca5(PO4)3⋅1OHNa"[n"H"], 0.0019, atol=0.0001)
    @test isapprox(mat"Ca5(PO4)3⋅1OHNa"[n"O"], 0.3959, atol=0.0001)
    
end