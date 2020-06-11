using Test

@testset "CustomMAC" begin
    @test NeXLCore.getcustommac(n"Si", n"O K-L3", :Henke1974) == 8790.0
    @test NeXLCore.getcustommac(n"U", n"B K-L2", :Henke1982) == 9020.0
    @test NeXLCore.getcustommac(n"Ta", n"N K-L3", :Bastin1988) == 15000.0

    cust1 = getcustommacs(:Bastin1988, true)
    cust2 = getcustommacs(:Bastin1988, false)

    @test !isnothing(findfirst(x->x==(n"Si", n"N K-L3", 17170.0), cust1))
    @test !isnothing(findfirst(x->x==(n"Si", n"N K-L2", 17170.0), cust1))
    @test isnothing(findfirst(x->x==(n"Si", n"N K-L2", 17170.0), cust2))
    @test !isnothing(findfirst(x->x==(n"Zr", n"N K-L3", 24000.0), cust1))

    clearusermacs()
    addcustommacs(:Bastin1988, true)
    
    @test mac(n"Si", n"N K-L3", UserMAC) == 17170.0
    @test mac(n"Zr", n"N K-L3", UserMAC) == 24000.0
    @test mac(n"Zr", n"O K-L3", UserMAC) == mac(n"Zr", n"O K-L3", FFASTDB)

    addcustommacs(:Bastin1989, true)
    @test mac(n"Zr", n"O K-L3", UserMAC) == 16200.0
    @test mac(n"Zr", n"O K-L2", UserMAC) == 16200.0
    @test mac(n"Pb", n"O K-L3", UserMAC) == 11000.0

    clearusermacs()
    @test mac(n"Zr", n"O K-L3", UserMAC) == mac(n"Zr", n"O K-L3", FFASTDB)
end
