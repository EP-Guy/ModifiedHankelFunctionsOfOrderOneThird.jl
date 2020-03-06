using Test
using ModifiedHankelFunctionsOfOrderOneThird

@testset "Zeros of h₂" begin
    h1, h2, h1p, h2p = modifiedhankel(complex(-1.16905371, 2.02486041))
    @test h1 ≈ complex(-0.34342750, 0.59483387)
    @test h1p ≈ complex(-0.06957489, -1.06099210)
    @test h2p ≈ complex(-1.83769223, 1.06099210)

    h1, h2, h1p, h2p = modifiedhankel(complex(-2.04397472, 3.54026807))
    @test h1 ≈ complex(0.29985265, -0.51936003)
    @test h1p ≈ complex(0.03621447, 1.21517637)
    @test h2p ≈ complex(2.10474721, -1.21517637)

    h1, h2, h1p, h2p = modifiedhankel(complex(-2.76027991, 4.78094505))
    @test h1 ≈ complex(-0.27833328, 0.48208739) atol=1e-7
    @test h1p ≈ complex(-0.02507370, -1.30912788)
    @test h2p ≈ complex(-2.26747600, 1.30912788)
end

@testset "Zeros of h₂′" begin
    h1, h2, h1p, h2p = modifiedhankel(complex(-0.50939649, 0.88230059))
    @test h1 ≈ complex(-0.63166599, -0.52691131)
    @test h2 ≈ complex(0, 1.62098891)
    @test h1p ≈ complex(0.89913968, 0)

    h1, h2, h1p, h2p = modifiedhankel(complex(-1.62409879, 2.81302162))
    @test h1 ≈ complex(0.53589957, 0.33980741)
    @test h2 ≈ complex(0, -1.26801270)
    @test h1p ≈ complex(-1.14943284, 0)

    h1, h2, h1p, h2p = modifiedhankel(complex(-2.41004961, 4.17432837))
    @test h1 ≈ complex(-0.49173751, -0.29946086) atol=1e-7
    @test h2 ≈ complex(0, 1.15117521)
    @test h1p ≈ complex(1.26609349, 0)
end

@testset "Spot checks" begin
    h1, h2, h1p, h2p = modifiedhankel(complex(2., 0.))
    @test h1 ≈ complex(0.60991262, 0.36822576)
    @test h1p ≈ complex(-0.59922801, 0.83306447)
    @test h2 ≈ complex(0.60991262, -0.36822576)
    @test h2p ≈ complex(-0.59922801, -0.83306447)

    h1, h2, h1p, h2p = modifiedhankel(complex(4.5, 2.2))
    @test h1 ≈ complex(-0.00170177, -0.00480825) atol=1e-7
    @test h1p ≈ complex(0.01152459, -0.00111977) atol=1e-7
    @test h2 ≈ complex(-7.26047008, 63.40667850)
    @test h2p ≈ complex(133.64890371, 44.66625637)

    h1, h2, h1p, h2p = modifiedhankel(complex(1.1, 4.1))
    @test h1 ≈ complex(-0.00202396, 0.00177285) atol=1e-7
    @test h1p ≈ complex(-0.00042308, -0.00566719) atol=1e-7
    @test h2 ≈ complex(-131.37854842, -8.34958396)
    @test h2p ≈ complex(-176.03846956, 196.31112124)
end

@testset "Asymptotic expansion" begin
    @test ModifiedHankelFunctionsOfOrderOneThird.C(1) ≈ 0.10416666666666667
    @test ModifiedHankelFunctionsOfOrderOneThird.C(2) ≈ 0.08355034722222222
    @test ModifiedHankelFunctionsOfOrderOneThird.C(5) ≈ 0.88162726744375765
    @test ModifiedHankelFunctionsOfOrderOneThird.C(7) ≈ 14.9957629868626
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(8)) ≈ 78.923013011587
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(11)) ≈ 24086.5496
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(14)) ≈ 17484377.
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(25)) ≈ 7.3900049415704853993e+19
end

@testset "LWPC spot checks" begin
    z = complex(24.6770630, 2.7361517)
    expz = complex(13.5990582,-80.0383987)

    # mh1lwpc = modhankelexpansion(z, false)
    h1, h2, h1p, h2p = ModifiedHankelFunctionsOfOrderOneThird.asymptotic(z)

    @test h1 ≈ complex(0.3822204,-0.0108716) *exp(-expz) rtol=1e-5
    @test h2 ≈ complex(0.3823441,-0.0102396) *exp(expz) rtol=1e-5
    @test h1p ≈ complex(-0.0548274,1.9051478) *exp(-expz) rtol=1e-5
    @test h2p ≈ complex(0.0503774,-1.9045298)*exp(expz) rtol=1e-5

    z = complex(21.7820282,2.7361517)
    expz = complex(12.7783260,-66.0632935)

    h1, h2, h1p, h2p = ModifiedHankelFunctionsOfOrderOneThird.asymptotic(z)

    @test h1 ≈ complex(0.3940974,-0.0127071)*exp(-expz) rtol=1e-5
    @test h2 ≈ complex(0.3942706,-0.0119274)*exp(expz) rtol=1e-5
    @test h1p ≈ complex(-0.0603090,1.8473313)*exp(-expz) rtol=1e-5
    @test h2p ≈ complex(0.0551328,-1.8465159)*exp(expz) rtol=1e-5

    z = complex(6.0799136,1.7321742)
    expz = complex(4.2853370,-8.3826580)

    h1, h2, h1p, h2p = ModifiedHankelFunctionsOfOrderOneThird.asymptotic(z)

    @test h1 ≈ complex(0.5353290,-0.0403038)*exp(-expz) rtol=1e-5
    @test h2 ≈ complex(0.5385904,-0.0340727)*exp(expz) rtol=1e-3
    @test h1p ≈ complex(-0.1057650,1.3544524)*exp(-expz) rtol=1e-5
    @test h2p ≈ complex(0.0823887,-1.3459417)*exp(expz) rtol=1e-3

    z = complex(74.5923615,-0.4041530)
    expz = complex(-3.4905479,-428.1735535)

    h1, h2, h1p, h2p = ModifiedHankelFunctionsOfOrderOneThird.asymptotic(z)

    @test h1 ≈ complex(0.2904782,3.4649431E-004)*exp(-expz) rtol=1e-3
    @test h2 ≈ complex(0.2904773,4.4042346E-004)*exp(expz) rtol=1e-3
    @test h1p ≈ complex(0.0028303,2.5087805)*exp(-expz) rtol=1e-3
    @test h2p ≈ complex(-0.0039661,-2.5087883)*exp(expz) rtol=1e-3
end
