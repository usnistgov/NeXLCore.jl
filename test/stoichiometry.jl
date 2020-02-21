using Test

using NeXLCore

@testset "Stoichiometry" begin

  # Check against DTSA-II
  @test isapprox(obystoichiometry(Dict(n"Mg"=>0.116567,n"Al"=>0.0490615,n"Si"=>0.211982, n"Ca"=>0.10899,n"Fe"=>0.0774196)), 0.4276, atol=0.0001)

  sio2 = asoxide(n"Si")
  @test name(sio2)=="SiO2"
  asio2 = atomicfraction("SiO2",Dict(n"Si"=>1,n"O"=>2))
  @test isapprox(sio2[n"Si"],asio2[n"Si"],atol=1.0e-6)
  @test isapprox(sio2[n"O"],asio2[n"O"],atol=1.0e-6)

  mix=asoxide(Dict(n"Si"=>0.4,n"Al"=>0.6))
  # Check against DTSA-II
  @test isapprox(mix[n"O"], 0.4955, atol=0.0001)
  @test isapprox(mix[n"Al"], 0.3176, atol=0.0001)
  @test isapprox(mix[n"Si"], 0.1870, atol=0.0001)

  amix = 0.4*atomicfraction("SiO2",Dict(n"Si"=>1,n"O"=>2))+0.6*atomicfraction("Al2O3",Dict(n"Al"=>2,n"O"=>3))
  @test isapprox(amix[n"O"], 0.4955, atol=0.0001)
  @test isapprox(amix[n"Al"], 0.3176, atol=0.0001)
  @test isapprox(amix[n"Si"], 0.1870, atol=0.0001)

  # Check against DTSA-II
  @test isapprox(amix[n"Si"], 0.1870,atol=0.0001)
  @test isapprox(amix[n"Al"], 0.3176,atol=0.0001)
  @test isapprox(amix[n"O"], 0.4955,atol=0.0001)
end
