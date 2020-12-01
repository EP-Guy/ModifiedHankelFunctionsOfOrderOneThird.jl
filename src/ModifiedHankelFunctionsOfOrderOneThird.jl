__precompile__(true)

module ModifiedHankelFunctionsOfOrderOneThird

using SpecialFunctions
using StaticArrays
using BetterExp
using Markdown
using Documenter

export modifiedhankel

# Number of terms to sum in both `powerseries` and `asymptotic` solutions
# Larger `NUMTERMS` generate successively smaller values, but can lead to floating point
# errors when multiplied against successively larger powers of `z`.
const NUMTERMS = 30

# Series for auxiliary functions `f`, `f′`, `g`, and `g′` in `powerseries`.
const a₀ = cbrt(2)/gamma(2/3)
a(m) = m == 0 ? a₀ : -a(m-1)/((3m-1)*3m)
const avec = @SVector [a(i) for i = 0:NUMTERMS]

const b₀ = cbrt(2)/(3^(2/3)*gamma(4/3))
b(m) = m == 0 ? b₀ : -b(m-1)/((3m+1)*3m)
const bvec = @SVector [b(i) for i = 0:NUMTERMS]

const c₀ = a₀/2
c(m) = m == 0 ? c₀ : -c(m-1)/((3m+2)*3m)
const cvec = @SVector [c(i) for i = 0:NUMTERMS]

const d₀ = b₀
d(m) = m == 0 ? d₀ : -d(m-1)/((3m-2)*3m)
const dvec = @SVector [d(i) for i = 0:NUMTERMS]

Cterm(m) = 9*(2m - 1)^2 - 4
C(m) = convert(Float64, prod(Cterm, 1:m)/(2^(4m)*3^m*factorial(m)))
# `big` required for `factorial`; also, this begins with 1 to simplify evalpoly below
const Cvec = @SVector [i == 0 ? 1.0 : C(big(i)) for i = 0:NUMTERMS]
# begin with 0 to simplify evalpoly; Cvec[i+1] to offset from the 0 term of Cvec
const Cveci = @SVector [i == 0 ? 0.0 : Cvec[i+1]*i for i = 0:NUMTERMS]

"""
    modifiedhankel(z)

Return ``h₁``, ``h₂``, ``h₁'``, and ``h₂'``, the first and second modified Hankel functions
of order 1/3 and their derivatives.

These functions solve Stokes' equation ``d²u/dz² + zu = 0`` as a power series in ``z`` for
`abs2(z) < 36` and an approximate asymptotic expansion otherwise. The asymptotic solution is
necessary because the ``z³ⁱ`` in the power series blows up as ``i → ∞``.

!!! warning

    For very large arguments `z`, the exponential function called within the asymptotic expansion
    may overflow and return `Inf`, `Nan`, or 0 values. No warnings will be thrown from this
    module or Base for finite precision issues.

For more information, see [^SCL1945].

See also: [`powerseries`](@ref), [`asymptotic`](@ref)

# Examples

```jldoctest
julia> h1, h2, h1prime, h2prime = modifiedhankel(complex(2.687, -0.648));
```

# References

[^SCL1945]:

    The Staff of the Computation Library (1945), *Tables of the modified Hankel function of
    order one-third and of their derivatives.* Cambridge, MA: Harvard University Press.
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

function checkinput(z::Complex{Float64})
    """
    asymptotic will return NaN or 0.0 because`exp` will overflow for Float64 arguments after
    being raised to the 3/2 power. In particular, the variable `e2` is:
        e2 = exp(2im/3*z^(3/2) - 5im*π/12)
    There is a similar expression for `e3`.

    Working backwards from `log(prevfloat(typemax(Float64))) ≈ 700`, we can solve for `z`
    and get `z = 51.653 - 89.47im` which has an `abs2` of approximately 10000.
    """

    arg = 2im/3*sqrt(z)^3
    testval = exp(arg)

    if isinf(testval) || isnan(testval) || (!iszero(arg) & iszero(testval))
        @warn "Special functions within asymptotic expansion will wrap or overflow on
            Float64 $z. Try using `big(z)`."
    end
end

"""
    powerseries(z)

Return ``h₁``, ``h₂``, ``h₁'``, and ``h₂'``, the first and second modified Hankel functions
of order 1/3 and their derivatives using a *power series* with 30 terms.

For more information, see [^SCL1945].

See also: [`modifiedhankel`](@ref), [`asymptotic`](@ref)

# References

[^SCL1945]:

    The Staff of the Computation Library (1945), *Tables of the modified Hankel function of
    order one-third and of their derivatives.* Cambridge, MA: Harvard University Press.
"""
function powerseries(z)
    z² = z^2
    z³ = z²*z  # z³

    # Auxiliary functions, beginning with zeroth terms (a₀, b₀, c₀, d₀)
    fval = evalpoly(z³, avec.data)
    gval = evalpoly(z³, bvec.data)
    f′val = evalpoly(z³, cvec.data)
    g′val = evalpoly(z³, dvec.data)

    gval *= z
    f′val *= -z²

    k1 = 1im*sqrt(3)/3*(gval - 2*fval)
    k2 = 1im*sqrt(3)/3*(g′val - 2*f′val)

    h1 = gval + k1
    h2 = gval - k1
    h1p = g′val + k2
    h2p = g′val - k2

    return h1, h2, h1p, h2p
end

"""
    asymptotic(z)

Return ``h₁``, ``h₂``, ``h₁'``, and ``h₂'``, the first and second modified Hankel functions
of order 1/3 and their derivatives using an *asymptotic expansion*, useful for `abs2(z) > 36`.

For more information, see [^SCL1945].

See also: [`modifiedhankel`](@ref), [`powerseries`](@ref)

# References

[^SCL1945]:

    The Staff of the Computation Library (1945), *Tables of the modified Hankel function of
    order one-third and of their derivatives.* Cambridge, MA: Harvard University Press.
"""
function asymptotic(z)
    α = cbrt(2)*3^(1/6)/sqrt(π)
    rootz = sqrt(z)
    rootz_cubed = rootz*z  # z^(3/2)

    zinv = inv(z)
    zterm = 1im/rootz_cubed  # im*z^(-3/2)
    negative_zterm = -zterm

    # First term manually specified
    s = evalpoly(negative_zterm, Cvec.data)
    t = evalpoly(zterm, Cvec.data)
    sp = evalpoly(negative_zterm, Cveci.data)
    tp = evalpoly(zterm, Cveci.data)

    k1 = -3/2*zterm*zinv  # -3/2*im*z^(-5/2)
    sp *= k1
    tp *= k1  # yes, same for both

    tmpa = α/sqrt(rootz)  # α*z^(-1/4)
    k2 = 2im/3*rootz_cubed
    tmpb = k2 - 5im*π/12
    tmpc = k2 + 11im*π/12
    k3 = zinv/4  # 1/4/z

    e2 = exp(tmpb)  # exp(-tmpb) = 1/exp(tmpb)
    e2inv = inv(e2)
    h1 = e2*s
    h2 = t*e2inv
    h1p = e2*(s*(1im*rootz - k3) + sp) # rootz*z^(-1/4) = z^(1/4)
    h2p = (t*(-1im*rootz - k3) + tp)*e2inv

    e3 = exp(tmpc)  # exp(-tmpc) = 1/exp(tmpc)
    e3inv = inv(e3)
    argz = angle(z)
    if -4π/3 < argz < 0
        h1 += t*e3inv
        h1p += (t*(-1im*rootz - k3) + tp)*e3inv
    else
        h2 += e3*s
        h2p += e3*(s*(1im*rootz - k3) + sp)
    end

    h1 *= tmpa
    h2 *= tmpa
    h1p *= tmpa
    h2p *= tmpa

    return h1, h2, h1p, h2p
end

end # module
