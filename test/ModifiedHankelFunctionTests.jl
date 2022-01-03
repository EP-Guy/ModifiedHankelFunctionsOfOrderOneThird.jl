using Test
using ModifiedHankelFunctionsOfOrderOneThird
using SpecialFunctions

# SCL1945:
# "Tables of the modified Hankel function of order one-third and of their derivatives"
# Harvard University Press, 1945.

@testset "Power series: zeros of h₂" begin
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

@testset "Power series: zeros of h₂′" begin
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

@testset "Power series: spot checks" begin
    h1, h2, h1p, h2p = modifiedhankel(complex(2.0, 0.0))
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

@testset "Comparison with Airy and Bessel functions" begin
    z = complex(1.345,-0.245)
    T = 2/3*z^(3/2)
    k = (3/2)^(2/3)*(1 - im*sqrt(3)/3)
    β = 2^(4/3)*3^(1/6)

    h1, h2, h1p, h2p = modifiedhankel(z)
    h1m, h2m, h1mp, h2mp = modifiedhankel(-z)
    h1c, h2c, h1cp, h2cp = modifiedhankel(conj(z))

    ai, bi = airyai(z), airybi(z)
    aim, bim = airyai(-z), airybi(-z)
    aip, bip = airyaiprime(z), airybiprime(z)

    # See SCL1945 (22), (25)
    @test h2 ≈ conj(h1c)
    @test h2p ≈ conj(h1cp)

    # SCL1945 (18), (19)
    @test h1 ≈ k*(aim - im*bim)
    @test h2 ≈ conj(k)*(aim + im*bim)

    # SCL1945 (20), (21)
    @test ai ≈ 1/2k*h1m + 1/(2*conj(k))*h2m
    @test bi ≈ im/2k*h1m - im/(2*conj(k))*h2m

    # SCL1945 (23), (24)
    @test aip ≈ -1/(2k)*h1mp - 1/(2*conj(k))*h2mp
    @test bip ≈ -im/(2k)*h1mp + im/(2*conj(k))*h2mp

    # SCL1945 (26)
    h1t = modifiedhankel(cis(-2π/3)*z)[1]
    h2t = modifiedhankel(cis(2π/3)*z)[2]
    @test -im*β*aim ≈ h1t
    @test -im*β*aim ≈ -h2t
    @test -im*β*aim ≈ -cis(π/3)*(cis(π/3)*h1 + h2)

    # SCL1945 (27)
    h1t = modifiedhankel(conj(z)*exp(2π*im/3))[1]
    @test h1t ≈ -h2c
    @test h1t ≈ -conj(h1)

    # SCL1945 (28)
    h1tp = modifiedhankel(z*exp(2π*im/3))[3]
    @test h1tp ≈ cis(π/3)*h2p
    @test h1tp ≈ cis(π/3)*conj(h1cp)

    # SCL1945
    @test besselj(1/3,T) ≈ 1/2*T^(-1/3)*(h1 + h2)  # (36)
    @test bessely(1/3,T) ≈ 1/2im*T^(-1/3)*(h1 - h2)  # (37)
    @test hankelh1(1/3,T) ≈ T^(-1/3)*h1  # (38)
    @test hankelh2(1/3,T) ≈ T^(-1/3)*h2  # (39)
    @test besselj(-1/3,T) ≈ 1/2*T^(-1/3)*(cis(π/3)*h1 + cis(-π/3)*h2)  # (40)
    @test bessely(-1/3,T) ≈ 1/2im*T^(-1/3)*(cis(π/3)*h1 - cis(-π/3)*h2)  # (41)
    @test hankelh1(-1/3,T) ≈ cis(π/3)*T^(-1/3)*h1  # (42)
    @test hankelh2(-1/3,T) ≈ cis(-π/3)*T^(-1/3)*h2  # (43)
    @test besselj(2/3,T) ≈ 1/2*(2/3)^(1/3)*T^(-2/3)*(cis(-2π/3)*h1p + cis(2π/3)*h2p)  # (44)
    @test bessely(2/3,T) ≈ 1/2im*(2/3)^(1/3)*T^(-2/3)*(cis(-2π/3)*h1p - cis(2π/3)*h2p)  # (45)
    @test hankelh1(2/3,T) ≈ (2/3)^(1/3)*cis(-2π/3)*T^(-2/3)*h1p  # (46)
    @test hankelh2(2/3,T) ≈ (2/3)^(1/3)*cis(2π/3)*T^(-2/3)*h2p  # (47)
    @test besselj(-2/3,T) ≈ 1/2*(2/3)^(1/3)*T^(-2/3)*(h1p + h2p)  # (48)
    @test bessely(-2/3,T) ≈ 1/2im*(2/3)^(1/3)*T^(-2/3)*(h1p - h2p)  # (49)
    @test hankelh1(-2/3,T) ≈ (2/3)^(1/3)*T^(-2/3)*h1p  # (50)
    @test hankelh2(-2/3,T) ≈ (2/3)^(1/3)*T^(-2/3)*h2p  # (51)
end

@testset "Asymptotic expansion" begin
    # From Table 2, SCL1945

    @test ModifiedHankelFunctionsOfOrderOneThird.C(1) == 0.10416666666666667
    @test ModifiedHankelFunctionsOfOrderOneThird.C(2) == 0.08355034722222222
    @test ModifiedHankelFunctionsOfOrderOneThird.C(3) == 0.12822657455632716
    @test ModifiedHankelFunctionsOfOrderOneThird.C(4) == 0.29184902646414046
    @test ModifiedHankelFunctionsOfOrderOneThird.C(5) == 0.88162726744375765

    @test ModifiedHankelFunctionsOfOrderOneThird.C(6) ≈ 3.321408281862768 atol=1e-14
    @test ModifiedHankelFunctionsOfOrderOneThird.C(7) ≈ 14.9957629868626 atol=1e-12
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(8)) ≈ 78.923013011587 atol=1e-11
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(9)) ≈ 474.451538868 atol=1e-8
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(10)) ≈ 3207.490091 atol=1e-5

    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(11)) ≈ 24086.5496 atol=1e-3
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(12)) ≈ 198923.12 atol=1e-1
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(13)) ≈ 1791902.0 atol=1
    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(14)) ≈ 17484377 atol=1

    @test ModifiedHankelFunctionsOfOrderOneThird.C(BigInt(25)) ≈ 7.3900049415704853993e+19
end
