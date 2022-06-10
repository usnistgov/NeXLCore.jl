using Test
using NeXLCore

@testset "CustomMAC" begin
    loadcustommac!(n"Si", n"O K-L3", "Henke1974")
    @test mac(n"Si", n"O K-L3") == 8790.0
    loadcustommac!(n"U", n"B K-L2", "Henke1982")
    @test mac(n"U", n"B K-L2") == 9020.0
    loadcustommac!(n"Ta", n"N K-L3", "Bastin1988")
    @test mac(n"Ta", n"N K-L3") == 15000.0
    @test_throws ErrorException loadcustommac!(n"W", n"F K-L3", "Bastin1988")
    
    loadcustommacs!("Bastin1988", 1:99)
    @test mac(n"Si", n"N K-L3") == 17170.0
    @test mac(n"Zr", n"N K-L3") == 24000.0
    # Not present in Bastin1988
    @test mac(n"Zr", n"O K-L3") == mac(n"Zr", energy(n"O K-L3"))

    loadcustommacs!("Bastin1989", 1:99)
    @test mac(n"Zr", n"O K-L3") == 16200.0
    @test mac(n"Pb", n"O K-L3") == 11000.0
    resetmac!(n"Ga", n"O K-L2")
    # Cleared...
    @test mac(n"Ga", n"O K-L2") == mac(n"Ga", energy(n"O K-L2"))

    resetmacs!()
    @test mac(n"Zr", n"O K-L3") == mac(n"Zr", energy(n"O K-L3"))
end
