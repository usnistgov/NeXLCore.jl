using Test

@testset "Element" begin
    @test configuration(n"Fe") == "1s² 2s² 2p⁶ 3s² 3p⁶ 4s² 3d⁶"
    @test element(20) == n"Ca"
    @test a(n"Ca")==40.0784
    @test all(zz->z(elements[zz])==zz,1:92)
    @test symbol(n"Pd")=="Pd"
    @test density(n"Pb")==11.34
    @test atomic_weight[n"Pb"].lo == 206.14
    @test atomic_weight[n"Pb"].hi == 207.94
    @test isapprox(atomic_weight[n"Ca"], UncertainValue(40.078, 0.004), atol=0.00001)
end