using Test

@testset "CustomMAC" begin
    @test mac(n"Si", n"O K-L3", :Henke1974) == 8790.0
    @test mac(n"U", n"B K-L2", :Henke1982) == 9020.0
    @test mac(n"Ta", n"N K-L3", :Bastin1988) == 15000.0
    @test ismissing(mac(n"W", n"F K-L3", :Bastin1988))

    cust1 = getcustommacs(:Bastin1988, true)
    cust2 = getcustommacs(:Bastin1988, false)

    @test !isnothing(findfirst(x -> x == (n"Si", n"N K-L3", 17170.0), cust1))
    @test !isnothing(findfirst(x -> x == (n"Si", n"N K-L2", 17170.0), cust1))
    @test isnothing(findfirst(x -> x == (n"Si", n"N K-L2", 17170.0), cust2))
    @test !isnothing(findfirst(x -> x == (n"Zr", n"N K-L3", 24000.0), cust1))

    clear_user_macs!()
    for ( el, cxr, mac) in getcustommacs(:Bastin1988, true)
        set_user_mac!(el, cxr, mac)
    end
    @test mac(n"Si", n"N K-L3") == 17170.0
    @test mac(n"Zr", n"N K-L3") == 24000.0
    @test mac(n"Zr", n"O K-L3") == mac(n"Zr", n"O K-L3", FFASTDB)

    for ( el, cxr, mac) in getcustommacs(:Bastin1989, true)
        set_user_mac!(el, cxr, mac)
    end
    @test mac(n"Zr", n"O K-L3") == 16200.0
    @test mac(n"Zr", n"O K-L2") == 16200.0
    @test mac(n"Pb", n"O K-L3") == 11000.0
    delete_user_mac!(n"Ga", n"O K-L2")
    @test mac(n"Ga", n"O K-L2") == mac(n"Ga", n"O K-L2", FFASTDB)

    clear_user_macs!()
    @test mac(n"Zr", n"O K-L3") == mac(n"Zr", n"O K-L3", FFASTDB)
end
