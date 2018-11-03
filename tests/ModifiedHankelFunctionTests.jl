using Test
using Traceur

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


Traceur.warnings(() -> ModifiedHankelFunctionsOfOrderOneThird.modifiedhankel(complex(1.1, 4.1)))


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




"""
Asymptotic expansion for ``abs(z) > 6``
"""
function modhankelexpansion(z, bothnegative::Bool)
    cap = [1.0416666666666666663e-01,  8.3550347222222222116e-02,
           1.2822657455632716019e-01,  2.9184902646414046315e-01,
           8.8162726744375764874e-01,  3.3214082818627675264e+00,
           1.4995762986862554546e+01,  7.8923013011586517530e+01,
           4.7445153886826431887e+02,  3.2074900908906619004e+03,
           1.7919020077753438063e+06,  1.7484377180034121023e+07,
           2.4086549640874004605e+04,  1.9892311916950979121e+05,
           1.8370737967633072978e+08,  2.0679040329451551508e+09,
           2.4827519375935888472e+10,  3.1669454981734887315e+11,
           4.2771126865134715582e+12,  6.0971132411392560749e+13,
           9.1486942234356396792e+14,  1.4413525170009350101e+16,
           2.3788844395175757942e+17,  4.1046081600946921885e+18,
           7.3900049415704853993e+19,  1.3859220004603943141e+21,
           2.7030825930275761623e+22,  5.4747478619645573335e+23,
           1.1498937014386333524e+25,  2.5014180692753603969e+26]

    @assert abs2(z) > 36

    # Asymptotic expansion
    α = complex(8.53667218838951e-1)
    zpower = complex(1)
    mpower = complex(1)
    sum1 = complex(1)
    sum2 = complex(1)
    sum3 = complex(0)
    sum4 = complex(0)
    rootz = sqrt(z)
    rootz_cubed = rootz*z
    zterm = im/rootz_cubed
    mterm = -zterm
    dm = complex(0)
    term3 = complex(1)

    last = false
    m = 1
    while !last & (m <= 30)
        zpower *= zterm
        mpower *= mterm
        dm += complex(1)
        term1 = cap[m]*zpower
        term2 = cap[m]*mpower

        # abs(term2/term3) >= 1 && (last = true)

        sum1 += term1
        sum2 += term2
        sum3 += term1*dm
        sum4 += term2*dm
        term4 = term2*dm

        # (abs(real(term4)) <= 1e-5abs(real(sum4))) &
        #     (abs(imag(term4)) <= 1e-5abs(imag(sum4))) && (last = true)
        #
        # if last
        #     println("LAST")
        # end

        term3 = term2
        m += 1
    end
    sum3 *= zterm*complex(-1.5)/z
    sum4 *= zterm*complex(-1.5)/z

    term1 = (complex(-0.25)-im*rootz_cubed)/z
    term2 = (complex(-0.25)+im*rootz_cubed)/z

    zh1 = sum2
    zh1p = sum2*term2 + sum4
    zh2 = sum1
    zh2p = sum1*term1 + sum3

    zexp = -im*2/3*rootz_cubed + im*π*5/12

    if real(z) < 0
        exp1 = exp(zexp)
        exp2 = exp(zexp-im*π*4/3)

        if imag(z) >= 0
            zh2 *= exp1
            zh2p *= exp1

            if !bothnegative
                zh2 += zh1/exp2
                zh2p += zh1p/exp2
            end

            zh1 /= exp1
            zh1p /= exp1
        else
            th2 = -zh1/exp2
            th2p = -zh1p/exp2

            zh1 = zh1/exp1 + zh2*exp2
            zh1p = zh1p/exp1 + zh2p*exp2

            if !bothnegative
                zh2 *= exp1
                zh2p *= exp1
            else
                zh2 = th2
                zh2p = th2p
            end
        end
        zexp = complex(0)
    end

    zterm = α/sqrt(rootz)
    zh1 *= zterm
    zh2 *= zterm
    zh1p *= zterm
    zh2p *= zterm

    return zh1, zh2, zh1p, zh2p, zexp
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
