__precompile__()

module ModifiedHankelFunctionsOfOrderOneThird

using SpecialFunctions
using StaticArrays
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
C(m)::Float64 = prod(Cterm, 1:m)/(2^(4m)*3^m*factorial(m))
const Cvec = @SVector [C(BigInt(i)) for i = 1:NUMTERMS]


"""
    modifiedhankel(z)

Return ``h₁``, ``h₂``, ``h₁'``, and ``h₂'``, the first and second modified Hankel functions
of order 1/3 and their derivatives.

These functions solve Stokes' equation ``d²u/dz² + zu = 0`` as a power series in ``z`` for
`abs2(z) < 36` and an approximate asymptotic expansion otherwise. The asymptotic solution is
necessary because the ``z³ⁱ`` in the power series blows up as ``i → ∞``.

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
    # Auxiliary functions
    # The zeroth terms (a₀, b₀, c₀, d₀)
    fval = oftype(z, avec[1])
    gval = oftype(z, bvec[1])
    f′val = oftype(z, cvec[1])
    g′val = oftype(z, dvec[1])

    z² = z^2
    zterm = z²*z # z³
    zpow = one(z)
    @inbounds for i = 2:NUMTERMS
        # TODO: Stop interations based on accuracy/convergence criteria
        zpow *= zterm  # for z^(3i)
        fval += avec[i]*zpow
        gval += bvec[i]*zpow
        f′val += cvec[i]*zpow
        g′val += dvec[i]*zpow
    end
    gval *= z
    f′val *= -z²

    isqrt3over3 = im*sqrt(3)/3
    k1 = isqrt3over3*(gval - 2*fval)
    k2 = isqrt3over3*(g′val - 2*f′val)

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
function asymptotic(z::T) where T
    α = cbrt(2)*3^(1/6)/sqrt(π)
    rootz = sqrt(z)
    rootz_cubed = rootz*z  # z^(3/2)

    zinv = inv(z)
    zterm = im/rootz_cubed  # im*z^(-3/2)
    negative_zterm = -zterm

    cT = complex(T)
    zpower = one(cT)
    negative_zpower = one(cT)
    s = one(cT)
    t = one(cT)
    sp = zero(cT)
    tp = zero(cT)
    @inbounds for i = 1:NUMTERMS
        # TODO: Stop interations based on accuracy/convergence criteria
        zpower *= zterm
        negative_zpower *= negative_zterm

        tmp1 = Cvec[i]*negative_zpower
        tmp2 = Cvec[i]*zpower

        s += tmp1
        t += tmp2
        sp += tmp1*i
        tp += tmp2*i
    end
    k1 = -3/2*zterm*zinv  # -3/2*im*z^(-5/2)
    sp *= k1
    tp *= k1  # yes, same for both

    tmpa = α/sqrt(rootz)  # α*z^(-1/4)
    k2 = 2/3*im*rootz_cubed
    tmpb = k2 - 5π*im/12
    tmpc = k2 + 11π*im/12
    k3 = zinv/4  # 1/4/z

    e2 = exp(tmpb)  # exp(-tmpb) = 1/exp(tmpb)
    e2inv = inv(e2)
    h1 = e2*s
    h2 = t*e2inv
    h1p = e2*(s*(im*rootz - k3) + sp) # rootz*z^(-1/4) = z^(1/4)
    h2p = (t*(-im*rootz - k3) + tp)*e2inv

    e3 = exp(tmpc)  # exp(-tmpc) = 1/exp(tmpc)
    e3inv = inv(e3)
    argz = angle(z)
    if -4π/3 < argz < 0
        h1 += t*e3inv
        h1p += (t*(-im*rootz - k3) + tp)*e3inv
    else
        h2 += e3*s
        h2p += e3*(s*(im*rootz - k3) + sp)
    end

    h1 *= tmpa
    h2 *= tmpa
    h1p *= tmpa
    h2p *= tmpa

    return h1, h2, h1p, h2p
end

end # module
