using Test
using NeXLCore

@testset "Iron" begin
    @test length(NeXLCore.eachedge(26)) == 10

    @test NeXLCore.hasedge(26, 1)
    @test NeXLCore.hasedge(26, 2)
    @test NeXLCore.hasedge(26, 3)
    @test NeXLCore.hasedge(26, 4)
    @test NeXLCore.hasedge(26, 5)
    @test NeXLCore.hasedge(26, 6)
    @test NeXLCore.hasedge(26, 7)
    @test NeXLCore.hasedge(26, 8)
    @test NeXLCore.hasedge(26, 9)
    @test !NeXLCore.hasedge(26, 11)

    @test all(i -> i in NeXLCore.eachedge(26), 1:10)
    @test !(11 in NeXLCore.eachedge(26))
    @test 2 in NeXLCore.eachedge(26)
    @test 7 in NeXLCore.eachedge(26)

    @test NeXLCore.edgeenergy(26, 1) == 7117.0
    @test NeXLCore.edgeenergy(26, 9) == 9.0
    @test NeXLCore.edgeenergy(26, 4) ≈ 713.0

    #@test isapprox(atomicweight(26), 55.845, atol = 0.003)

    #@test density(26) ≈ 7.860

    #@test EV(26) ≈ 7.53493E+05

    #@test crosssectionfactor(26) ≈ 9.27362E+01 * 1.0e-24

    #@test relativisticcorrection(26)[1] ≈ -1.0812E-01
    #@test relativisticcorrection(26)[2] ≈ -6.7800E-02

    #@test nuclearthompsoncorrection(26) ≈ -6.6403E-03

    @test isapprox(mac(26, 30.0), 106592.1003, rtol = 0.001)
    @test isapprox(mac(26, 102.0), 46393.2381, rtol = 0.001)
    @test isapprox(mac(26, 304.3), 10243.4145, rtol = 0.001)
    @test isapprox(mac(26, 1234.1), 5297.2337, rtol = 0.001)
    @test isapprox(mac(26, 2342.0), 1013.5935, rtol = 0.001)
    @test isapprox(mac(26, 4523.0), 167.5849, rtol = 0.001)
    @test isapprox(mac(26, 12340.0), 95.1269, rtol = 0.001)
    @test isapprox(mac(26, 34000.0), 5.2825, rtol = 0.001)
    @test isapprox(mac(26, 9630.0), 186.2521, rtol = 0.001)
    @test isapprox(mac(26, 201000.0), 0.0233, rtol = 0.001)
end


@testset "Oxygen" begin
    @test length(NeXLCore.eachedge(8)) == 4

    @test NeXLCore.hasedge(8, 1)
    @test NeXLCore.hasedge(8, 2)
    @test NeXLCore.hasedge(8, 3)
    @test NeXLCore.hasedge(8, 4)
    @test !NeXLCore.hasedge(8, 5)
    @test !NeXLCore.hasedge(8, 6)
    @test !NeXLCore.hasedge(8, 7)
    @test !NeXLCore.hasedge(8, 8)
    @test !NeXLCore.hasedge(8, 9)
    @test !NeXLCore.hasedge(8, 10)

    @test NeXLCore.edgeenergy(8, 1) == 5.38000E-01 * 1000.0
    @test NeXLCore.edgeenergy(8, 2) ≈ 2.848000E-02 * 1000.0
    @test NeXLCore.edgeenergy(8, 4) ≈ 13.620000E-03 * 1000.0

    #@test isapprox(atomicweight(8), 15.99940, atol = 0.003)

    #@test density(8) ≈ 1.3310E-03

    #@test EV(8) ≈ 2.63012E+06

    #@test crosssectionfactor(8) ≈ 2.65676e1 * 1.0e-24

    #@test relativisticcorrection(8)[1] ≈ -7.7133E-03
    #@test relativisticcorrection(8)[2] ≈ -4.2000E-03

    #@test nuclearthompsoncorrection(8) ≈ -2.1944E-03

    @test isapprox(mac(8, 30.0), 303584.6463, rtol = 0.001)

    @test isapprox(mac(8, 30.0), 303584.6463, rtol = 0.001)
    @test isapprox(mac(8, 102.0), 49996.7695, rtol = 0.001)
    @test isapprox(mac(8, 304.3), 4191.9561, rtol = 0.001)
    @test isapprox(mac(8, 1234.1), 2491.7539, rtol = 0.001)
    @test isapprox(mac(8, 2342.0), 421.9132, rtol = 0.001)
    @test isapprox(mac(8, 4523.0), 59.6995, rtol = 0.001)
    @test isapprox(mac(8, 12340.0), 2.7359, rtol = 0.001)
    @test isapprox(mac(8, 34000.0), 0.1055, rtol = 0.001)
    @test isapprox(mac(8, 9630.0), 6.0225, rtol = 0.001)
    @test isapprox(mac(8, 201000.0), 0.00030926, rtol = 0.001)
end

@testset "Lead" begin
    @test length(NeXLCore.eachedge(82)) == 24

    @test NeXLCore.hasedge(82, 1)
    @test NeXLCore.hasedge(82, 2)
    @test NeXLCore.hasedge(82, 3)
    @test NeXLCore.hasedge(82, 4)
    @test NeXLCore.hasedge(82, 5)
    @test NeXLCore.hasedge(82, 6)
    @test NeXLCore.hasedge(82, 7)
    @test NeXLCore.hasedge(82, 8)
    @test !NeXLCore.hasedge(82, 23)
    @test !NeXLCore.hasedge(82, 25)
    @test NeXLCore.hasedge(82, 26)
    @test NeXLCore.hasedge(82, 27)

    @test isapprox(mac(82, 30.0), 74494.9013, rtol = 0.001)
    @test isapprox(mac(82, 102.0), 11406.5505, rtol = 0.001)
    @test isapprox(mac(82, 304.3), 12149.8231, rtol = 0.001)
    @test isapprox(mac(82, 1234.1), 3192.3245, rtol = 0.001)
    @test isapprox(mac(82, 2342.0), 841.2044, rtol = 0.001)
    @test isapprox(mac(82, 4523.0), 893.4338, rtol = 0.001)
    @test isapprox(mac(82, 12340.0), 69.5927, rtol = 0.001)
    @test isapprox(mac(82, 34000.0), 20.0928, rtol = 0.001)
    @test isapprox(mac(82, 9630.0), 134.1960, rtol = 0.001)
    @test isapprox(mac(82, 201000.0), 0.8207, rtol = 0.001)
end

@testset "MACpe" begin # generated using NIST DTSA-II
    @test isapprox(mac(1, 1825.6371), 0.9074, rtol = 0.001)
    @test isapprox(mac(2, 67143.2841), 2.768340e-5, rtol = 0.001)
    @test isapprox(mac(3, 52.8894), 6341.1206, rtol = 0.001)
    @test isapprox(mac(4, 323.5202), 13139.1498, rtol = 0.001)
    @test isapprox(mac(5, 72.4552), 22406.1014, rtol = 0.001)
    @test isapprox(mac(6, 131.0993), 12738.7273, rtol = 0.001)
    @test isapprox(mac(7, 54.3167), 111770.5028, rtol = 0.001)
    @test isapprox(mac(8, 149.9255), 22397.5204, rtol = 0.001)
    @test isapprox(mac(9, 34281.3188), 0.1490, rtol = 0.001)
    @test isapprox(mac(10, 202.6272), 23577.7535, rtol = 0.001)
    @test isapprox(mac(11, 91714.0014), 0.0129137, rtol = 0.001)
    @test isapprox(mac(12, 51.6654), 57373.9013, rtol = 0.001)
    @test isapprox(mac(13, 93123.8498), 0.0230, rtol = 0.001)
    @test isapprox(mac(14, 80174.5908), 0.0503, rtol = 0.001)
    @test isapprox(mac(15, 19484.3456), 5.1816, rtol = 0.001)
    @test isapprox(mac(16, 51.6164), 12953.7948, rtol = 0.001)
    @test isapprox(mac(17, 92.6523), 14039.6832, rtol = 0.001)
    @test isapprox(mac(18, 286446.7751), 0.0019365, rtol = 0.001)
    @test isapprox(mac(19, 1282.5060), 2023.8581, rtol = 0.001)
    @test isapprox(mac(20, 58.1960), 21012.6801, rtol = 0.001)
    @test isapprox(mac(21, 2663.9168), 380.4134, rtol = 0.001)
    @test isapprox(mac(22, 52.0073), 12148.6626, rtol = 0.001)
    @test isapprox(mac(23, 51.4227), 56667.2105, rtol = 0.001)
    @test isapprox(mac(24, 22456.4234), 13.8940, rtol = 0.001)
    @test isapprox(mac(25, 11995.0019), 91.0401, rtol = 0.001)
    @test isapprox(mac(26, 55.5593), 106613.1135, rtol = 0.001)
    @test isapprox(mac(27, 22266.9550), 19.8706, rtol = 0.001)
    @test isapprox(mac(28, 205.8346), 22273.6361, rtol = 0.001)
    @test isapprox(mac(29, 1039.2688), 9356.0732, rtol = 0.001)
    @test isapprox(mac(30, 485.3171), 6755.4381, rtol = 0.001)
    @test isapprox(mac(31, 56.2430), 75027.5850, rtol = 0.001)
    @test isapprox(mac(32, 164.5878), 47424.7131, rtol = 0.001)
    @test isapprox(mac(33, 14923.1883), 97.4613, rtol = 0.001)
    @test isapprox(mac(34, 442.2493), 12244.0970, rtol = 0.001)
    @test isapprox(mac(35, 30402.3441), 16.1595, rtol = 0.001)
    @test isapprox(mac(36, 1763.8442), 4211.7645, rtol = 0.001)
    @test isapprox(mac(37, 446.8578), 16086.7643, rtol = 0.001)
    @test isapprox(mac(38, 96545.4058), 0.7178, rtol = 0.001)
    @test isapprox(mac(39, 822.5155), 5704.5987, rtol = 0.001)
    @test isapprox(mac(40, 52279.3424), 5.0292, rtol = 0.001)
    @test isapprox(mac(41, 127578.1715), 0.4067, rtol = 0.001)
    @test isapprox(mac(42, 18169.7086), 15.0481, rtol = 0.001)
    @test isapprox(mac(43, 119.6601), 5432.8402, rtol = 0.001)
    @test isapprox(mac(44, 3314.3650), 1661.1260, rtol = 0.001)
    @test isapprox(mac(45, 613.7420), 15720.1429, rtol = 0.001)
    @test isapprox(mac(46, 297.1361), 2936.5826, rtol = 0.001)
    @test isapprox(mac(47, 1798.6100), 1682.0306, rtol = 0.001)
    @test isapprox(mac(48, 201700.6027), 0.1723, rtol = 0.001)
    @test isapprox(mac(49, 68.2226), 74163.3914, rtol = 0.001)
    @test isapprox(mac(50, 227.3590), 6691.7553, rtol = 0.001)
    @test isapprox(mac(51, 91067.8017), 1.9996, rtol = 0.001)
    @test isapprox(mac(52, 7352.2301), 327.2407, rtol = 0.001)
    @test isapprox(mac(53, 271.7124), 4697.0153, rtol = 0.001)
    @test isapprox(mac(54, 96.1322), 31927.3174, rtol = 0.001)
    @test isapprox(mac(55, 452.1765), 4740.8624, rtol = 0.001)
    @test isapprox(mac(56, 2334.6617), 1522.5183, rtol = 0.001)
    @test isapprox(mac(57, 103.9157), 39300.3286, rtol = 0.001)
    @test isapprox(mac(58, 4223.4834), 392.9113, rtol = 0.001)
    @test isapprox(mac(59, 140717.3325), 0.8885, rtol = 0.001)
    @test isapprox(mac(60, 261513.7888), 0.1654, rtol = 0.001)
    @test isapprox(mac(61, 156.0607), 11951.7922, rtol = 0.001)
    @test isapprox(mac(62, 106138.7999), 2.2293, rtol = 0.001)
    @test isapprox(mac(63, 189663.2918), 0.4635, rtol = 0.001)
    @test isapprox(mac(64, 6110.7096), 204.0632, rtol = 0.001)
    @test isapprox(mac(65, 44195.9461), 4.8034, rtol = 0.001)
    @test isapprox(mac(66, 20052.2771), 44.3789, rtol = 0.001)
    @test isapprox(mac(67, 71.5978), 28290.5339, rtol = 0.001)
    @test isapprox(mac(68, 3221.8084), 1245.9590, rtol = 0.001)
    @test isapprox(mac(69, 58.3034), 31253.1709, rtol = 0.001)
    @test isapprox(mac(70, 84.7028), 26314.9599, rtol = 0.001)
    @test isapprox(mac(71, 95.3856), 25954.8160, rtol = 0.001)
    @test isapprox(mac(72, 10678.1649), 184.5189, rtol = 0.001)
    @test isapprox(mac(73, 86298.9086), 5.8990, rtol = 0.001)
    @test isapprox(mac(74, 59.5091), 19804.6389, rtol = 0.001)
    @test isapprox(mac(75, 1742.5269), 1093.5987, rtol = 0.001)
    @test isapprox(mac(76, 874.9993), 4872.1773, rtol = 0.001)
    @test isapprox(mac(77, 778.1277), 6112.9332, rtol = 0.001)
    @test isapprox(mac(78, 77.9444), 15733.9354, rtol = 0.001)
    @test isapprox(mac(79, 53.5598), 19230.5901, rtol = 0.001)
    @test isapprox(mac(80, 313.9529), 15487.9170, rtol = 0.001)
    @test isapprox(mac(81, 148.5939), 4791.9895, rtol = 0.001)
    @test isapprox(mac(82, 169292.0552), 1.2947, rtol = 0.001)
    @test isapprox(mac(83, 158.7399), 4193.2797, rtol = 0.001)
    @test isapprox(mac(84, 11247.1253), 97.4316, rtol = 0.001)
    @test isapprox(mac(85, 64952.8990), 4.0008, rtol = 0.001)
    @test isapprox(mac(86, 24619.4642), 53.9372, rtol = 0.001)
    @test isapprox(mac(87, 52.2036), 4222.7289, rtol = 0.001)
    @test isapprox(mac(88, 51.5550), 5607.8556, rtol = 0.001)
    @test isapprox(mac(89, 1035.3394), 5796.1285, rtol = 0.001)
    @test isapprox(mac(90, 142.6819), 6556.7729, rtol = 0.001)
    @test isapprox(mac(91, 51.2632), 8270.9798, rtol = 0.001)
    @test isapprox(mac(92, 51.8461), 13886.1866, rtol = 0.001)
end