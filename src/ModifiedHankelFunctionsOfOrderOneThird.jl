__precompile__()

module ModifiedHankelFunctionsOfOrderOneThird

using SpecialFunctions
using StaticArrays
using Markdown

export modifiedhankel

# Number of terms to sum in both `powerseries` and `asymptotic` solutions
# Larger `numterms` generate successively smaller values, but can lead to floating point
# errors when multiplied against successively larger powers of `z`.
const numterms = 30

# Series for auxiliary functions `f`, `f′`, `g`, and `g′` in `powerseries`.
const a₀ = cbrt(2)/gamma(2/3)
a(m) = m == 0 ? a₀ : -a(m-1)/((3m-1)*3m)
const avec = @SVector [a(i) for i = 0:numterms]

const b₀ = cbrt(2)/(3^(2/3)*gamma(4/3))
b(m) = m == 0 ? b₀ : -b(m-1)/((3m+1)*3m)
const bvec = @SVector [b(i) for i = 0:numterms]

const c₀ = a₀/2
c(m) = m == 0 ? c₀ : -c(m-1)/((3m+2)*3m)
const cvec = @SVector [c(i) for i = 0:numterms]

const d₀ = b₀
d(m) = m == 0 ? d₀ : -d(m-1)/((3m-2)*3m)
const dvec = @SVector [d(i) for i = 0:numterms]

Cterm(m) = 9(2m - 1)^2 - 4
C(m)::Float64 = prod(Cterm, 1:m)/(2^(4m)*3^m*factorial(m))
const Cvec = @SVector [C(BigInt(i)) for i = 1:numterms]


"""
Return h₁, h₂, h₁′, and h₂′, the first and second modified Hankel functions of order 1/3 and
their derivatives.

These functions solve Stokes' equation ``d²u/dz² + zu = 0`` as a power series in ``z`` for
`abs(z) < 6` and an approximate asymptotic expansion otherwise. The asymptotic solution is
necessary because the ``z³ⁱ`` in the power series blows up as ``i → ∞``.

Based on Tables of the Modified Hankel Functions...

See also: [powerseries](@ref), [asymptotic](@ref)
"""
function modifiedhankel(z)
    if abs2(z) < 36
        # Power series solution
        h1, h2, h1p, h2p = powerseries(z)
    else
        # Asymptotic expansion
        h1, h2, h1p, h2p = asymptotic(z)
    end

    return h1, h2, h1p, h2p
end

"""
    powerseries(z)

Return h₁, h₂, h₁′, and h₂′, the first and second modified Hankel functions of order 1/3 and
their derivatives using a *power series*.

These functions solve Stokes' equation ``d²u/dz² + zu = 0`` as a power series in ``z``,
valid in the entire complex plane. Based on Tables of the Modified Hankel Functions...

See also: [modifiedhankel](@ref), [asymptotic](@ref)
"""
function powerseries(z)
    # Auxiliary functions
    fval = zero(z)
    gval = zero(z)
    f′val = zero(z)
    g′val = zero(z)

    # The zeroth terms (a₀, b₀, c₀, d₀)
    fval += avec[1]
    gval += bvec[1]
    f′val += cvec[1]
    g′val += dvec[1]

    zterm = z^3
    zpow = one(z)
    @inbounds for i = 1:numterms
        zpow *= zterm  # for z^(3i)
        fval += avec[i+1]*zpow
        gval += bvec[i+1]*zpow
        f′val += cvec[i+1]*zpow
        g′val += dvec[i+1]*zpow
    end
    gval *= z
    f′val *= -z^2

    isqrt3over3 = im*sqrt(3)/3
    h1 = gval + isqrt3over3*(gval - 2fval)
    h2 = gval - isqrt3over3*(gval - 2fval)
    h1p = g′val + isqrt3over3*(g′val - 2f′val)
    h2p = g′val - isqrt3over3*(g′val - 2f′val)

    return h1, h2, h1p, h2p
end


doc"""
    asymptotic(z)

Return h₁, h₂, h₁′, and h₂′, the first and second modified Hankel functions of order 1/3 and
their derivatives using an *asymptotic expansion*.

Two separate functions, valid across separate ranges of ``arg(z)`` are necessary to cover
the full region of the complex plane. This is an example of Stokes' phenomenon, the existence
of two expressions of different forms which represent asymptotically the same integral function.

Based on Modified Hankel Functions ...

For ``h₁`` on ``-2π/3 < \arg z < 4π/3``:
```math
h₁(z) ≈ α z^{-1/4} e^{2/3 i z^{3/2} - 5πi/12} \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right)
h₁′(z) ≈ α i z^{1/4} e^{2/3 i z^{3/2} - 5πi/12} \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) -
    α/4 z^{-5/4} e^{2/3 i z^{3/2} - 5πi/12} \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) -
    3/2 α z^{-1/4} e^{2/3 i z^{3/2} - 5πi/12} \left( \sum_{m=1} (-i)^m m C_m z^{-3m/2 - 1} \right)
```
where
```math
C_m = \frac{(9-4)(81-4)\cdots (9[2m-1]^2-4)}{2^{4m}3^m m!}
```
and
```math
h₁(z) ≈ h₁(z) [above] + α z^{-1/4} e^{-2/3 i z^{3/2} - 11πi/12} \left( 1 + \sum_{m_1} (i)^m C_m z^{-3m/2} \right)
h₁′(z) ≈ h₁′(z) [above] - α z^{1/4} e^{-2/3 i z^{3/2} - 11πi/12} \left( 1 + \sum_{m_1} (i)^m C_m z^{-3m/2} \right) -
    α/4 z^{-5/4} e^{-2/3 i z^{3/2} - 11πi/12} \left( 1 + \sum_{m_1} (i)^m C_m z^{-3m/2} \right) -
    3/2 α z^{-1/4} e^{-2/3 i z^{3/2} - 11πi/12} \left( \sum_{m=1} i^m m C_m z^{-3m/2 - 1} \right)
```
for ``-4π/3 < \arg z < 0``.

And for ``h₂`` on ``-4π/3 < \arg z < 2π/3``:
```math
h₂(z) ≈ α z^{-1/4} e^{-2/3 i z^{3/2} + 5πi/12} \left( 1 + \sum_{m=1} (i)^m C_m z^{-3m/2} \right)
h₂′(z) ≈ -α z^{1/4} e^{-2/3 i z^{3/2} + 5πi/12} \left( 1 + \sum_{m=1} (i)^m C_m z^{-3m/2} \right) -
    α/4 z^{-5/4} e^{-2/3 i z^{3/2} + 5πi/12} \left( 1 + \sum_{m=1} (i)^m C_m z^{-3m/2} \right) -
    3/2 α z^{-1/4} e^{-2/3 i z^{3/2} + 5πi/12} \left( i^m m C_m z^{-3m/2 - 1} \right)
```
and
```math
h₂(z) ≈ h₂(z) [above] + α z^{-1/4} e^{2/3 i z^{3/2} + 11πi/12} \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right)
h₂′(z) ≈ h₂′(z) [above] + α z^{1/4} e^{2/3 i z^{3/2} + 11πi/12} \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) -
    α/4 z^{-5/4} e^{2/3 i z^{3/2} + 11πi/12} \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) -
    3/2 α z^{-1/4} e^{2/3 i z^{3/2} + 11πi/12} \left( \sum_{m=1} (-i)^m m C_m z^{-3m/2 - 1} \right)
```
for ``0 < \arg z < 4π/3``.

See also: [modifiedhankel](@ref), [powerseries](@ref)
"""
function asymptotic(z)
    α = cbrt(2)*3^(1/6)/sqrt(π)
    rootz = sqrt(z)
    rootz_cubed = rootz*z  # z^(3/2)

    zterm = im/rootz_cubed  # im*z^(-3/2)
    negative_zterm = -zterm

    zpower = one(z)
    negative_zpower = one(z)
    s = one(z)
    t = one(z)
    sp = zero(z)
    tp = zero(z)
    @inbounds for i = 1:numterms
        zpower *= zterm
        negative_zpower *= negative_zterm

        tmp1 = Cvec[i]*negative_zpower
        tmp2 = Cvec[i]*zpower

        s += tmp1
        t += tmp2
        sp += tmp1*i
        tp += tmp2*i
    end
    sp *= -3/2*zterm/z  # -3/2*im*z^(-5/2)
    tp *= -3/2*zterm/z  # yes, same for both

    tmp1 = α/sqrt(rootz)  # α*z^(-1/4)
    tmp2 = 2/3*im*rootz_cubed - 5π*im/12
    tmp3 = 2/3*im*rootz_cubed + 11π*im/12

    h1 = exp(tmp2)*s
    h2 = exp(-tmp2)*t
    h1p = exp(tmp2)*(s*(im*rootz - 1/4/z) + sp) # rootz*z^(-1/4) = z^(1/4)
    h2p = exp(-tmp2)*(t*(-im*rootz - 1/4/z) + tp)

    argz = angle(z)
    if -4π/3 < argz < 0
        h1 += exp(-tmp3)*t
        h1p += exp(-tmp3)*(t*(-im*rootz - 1/4/z) + tp)
    else
        h2 += exp(tmp3)*s
        h2p += exp(tmp3)*(s*(im*rootz - 1/4/z) + sp)
    end

    h1 *= tmp1
    h2 *= tmp1
    h1p *= tmp1
    h2p *= tmp1

    return h1, h2, h1p, h2p
end

end # module
