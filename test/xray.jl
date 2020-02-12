using Test
using NeXLCore

@testset "X-ray" begin

        e1, e2, e3 = element(20), element(84), element(92)

        @testset "Elements" begin

                @test z(e1)==20
                @test z(e2)==84
                @test z(e3)==92

                @test isapprox(a(e1), 40.08,atol=0.01)
                @test isapprox(a(e2),209.00,atol=0.01)
                @test isapprox(a(e3),238.03,atol=0.01)

                @test symbol(e1)=="Ca"
                @test symbol(e2)=="Po"
                @test symbol(e3)=="U"

                @test name(e1)=="Calcium"
                @test name(e2)=="Polonium"
                @test name(e3)=="Uranium"

                @test e1 ≠ e2
                @test e1 == element(20)
                @test e2 == element(84)
                @test e3 == element(92)

                @test e1 < e2
                @test e1 < e3
                @test e2 < e3

                @test e2 >= e1
                @test e2 >= e2
                @test e3 >= e2
                @test e3 >= e1

                @test e2 > e1
                @test e3 > e2
                @test e3 > e1
        end

        @testset "Shells" begin
                @test length(ksubshells)==1
                @test length(lsubshells)==3
                @test length(msubshells)==5
                @test length(nsubshells)==7
                @test length(osubshells)==9

                @test ksubshells == map(SubShell, ( "K", ))
                @test lsubshells == map(SubShell, ( "L1", "L2", "L3" ))
                @test msubshells == map(SubShell, ( "M1", "M2", "M3", "M4", "M5" ))
                @test nsubshells == map(SubShell, ( "N1", "N2", "N3", "N4", "N5", "N6", "N7" ))
                @test osubshells == map(SubShell, ( "O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9", ))

                @test shell(SubShell("K"))=='K'
                @test shell(SubShell("L2"))=='L'
                @test shell(SubShell("M4"))=='M'
                @test shell(SubShell("N7"))=='N'
                @test shell(SubShell("O1"))=='O'

                @test capacity(n"K1")==2
                @test capacity(n"L2")==2
                @test capacity(n"M5")==6
                @test capacity(n"N6")==6
                @test capacity(n"O9")==10

                as1, as2, as3 = AtomicSubShell(e1.number, subshell("L3")), AtomicSubShell(e2.number, subshell("M5")), AtomicSubShell(e3.number, subshell("K"))

                @test repr(as1)=="Ca L3"
                @test repr(as2)=="Po M5"
                @test repr(as3)=="U K"

                @test as1==n"Ca L3"
                @test as2==n"Po M5"
                @test as3==n"U K"

                @test shell(as1)=='L'
                @test shell(as2)=='M'
                @test shell(as3)=='K'

                @test n"Ca K" < n"Ca L1"
                @test n"Ca K" < n"Ti K"
                @test as1 < as2
                @test !(as2 < as1)
                @test !(as1 < as1)

                @test n"Ca K" == n"Ca K"
                @test !(n"Ca K" == n"Ca L1")
        end

        @testset "Transitions" begin
                @test n"K-L3" in alltransitions
                @test n"L3-M5" in alltransitions
                #@test_throws ArgumentError n"M5-L3"
                #@test_throws ArgumentError n"K-N7"

                @test shell(n"K-L3")=='K'
                @test shell(n"L3-M5")=='L'
                @test shell(n"M5-N7")=='M'

                @test all(tr->shell(tr)=='K',ktransitions)
                @test all(tr->shell(tr)=='L',ltransitions)
                @test all(tr->shell(tr)=='M',mtransitions)
                @test all(tr->shell(tr)=='N',ntransitions)

                @test length(alltransitions)==sum(length,(ktransitions,ltransitions,mtransitions,ntransitions, otransitions))

                @test n"K-L3"==n"K-L3"
                @test n"K-L3" ≠ n"K-L2"

                @test transition(n"Fe L3-M5")==n"L3-M5"
                @test element(n"Fe L3-M5")==n"Fe"

                @test isapprox(energy(n"Fe",n"L3-M5"),705.0,atol=1.0)
                @test isapprox(energy(n"U",n"L3-M5"),13614.6,atol=1.0)
                @test energy(n"Fe",n"L3-M5")==energy(n"Fe L3-M5")
                @test isapprox(energy(n"C K-L2"),277.4,atol=1.0)

                # from https://www.physics.nist.gov/cgi-bin/XrayTrans/search.pl?element=Fe&trans=KL3&lower=&upper=&units=A
                @test isapprox(λ(n"Fe K-L3"), 1.936e-8, atol = 1e-11) #
                @test isapprox(λ(n"Cu L3-M5"), 13.336e-8, atol=1e-11

                @test isapprox(energy(n"Fe L3"),708.0,atol=1.0)
                @test energy(n"Fe L3")==energy(n"Fe",n"L3")

                @test has(n"C",n"K-L2")
                @test !has(n"Li",n"K-L3")
                @test !has(n"C",n"L3-M5")
                @test has(n"U",n"L3-M5")
                @test has(n"U",n"M5-N7")

                @test all(tr->shell(tr)=='L',characteristic(n"Fe",ltransitions))
                @test all(tr->shell(tr)=='K',characteristic(n"Fe",ktransitions))

                @test length(characteristic(n"Fe",ltransitions,0.0))==12
                @test length(characteristic(n"Fe",ltransitions,0.1))==3
                @test length(characteristic(n"Fe",ltransitions,0.01))==7
        end
end
