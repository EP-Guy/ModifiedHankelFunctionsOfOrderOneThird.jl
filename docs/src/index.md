# Modified Hankel Functions of Order One Third

_Calculate modified Hankel functions of order one-third and their derivatives._

These special functions are solutions to Stokes' differential equation:
```math
\frac{\mathrm{d}^2u}{\mathrm{d}z^2} + zu = 0
```

The modified Hankel functions are linearly related to Airy functions but are not
as well known. One application of the modified Hankel functions of order one-third
is as a solution to the wave equation for electromagnetic fields between
the boundaries of the [earth-ionosphere waveguide](https://en.wikipedia.org/wiki/Earth%E2%80%93ionosphere_waveguide). It may also occur in other instances of
diffraction and refraction of waves.

The direct computation of solutions to Stokes' equation is preferred to
using Airy or Hankel functions. Unlike Bessel's equation, Stokes' equation has
no singularity in the finite complex plane and its solutions are single-valued
[^SCL1945].

## Usage

```julia
using ModifiedHankelFunctionsOfOrderOneThird

h1, h2, h1prime, h2prime = modifiedhankel(z)
```

## Solutions

Two solution approaches are used, as presented in [^SCL1945]. If `abs2(z) < 36`,
a power series solution is used. Otherwise, an asymptotic expansion is performed
because of floating point limits in the power series.

### Power Series Solution

``h₁``, ``h₂``, ``h₁'``, and ``h₂'`` are computed from auxiliary functions
```math
h_1(z) = g + \frac{i\sqrt{3}}{3}(g - 2f), \qquad h_2(z) = g - \frac{i\sqrt{3}}{3}(g - 2f)
```
where
```math
f(z) = a_0 + a_1z^3 + a_2z^6 + \cdots + a_mz^{3m} + \cdots
```
```math
g(z) = z(b_0 + b_1z^3 + b_2z^6 + \cdots + b_mz^{3m} + \cdots)
```
```math
f'(z) = -z^2(c_0 + c_1z^3 + c_2z^6 + \cdots + c_mz^{3m} + \cdots)
```
```math
g'(z) = d_0 + d_1z^3 + d_2z^6 + \cdots + d_mz^{3m} + \cdots
```
where
```math
a_0 = \frac{2^{1/3}}{\Gamma\left(\frac{2}{3}\right)} \qquad a_m = -\frac{a_{m-1}}{(3m-1)3m}
```
```math
b_0 = \frac{2^{1/3}}{3^{2/3}\Gamma\left(\frac{4}{3}\right)} \qquad b_m = -\frac{b_{m-1}}{3m(3m+1)}
```
```math
c_0 = \frac{a_0}{2} \qquad c_m = - \frac{c_{m-1}}{3m(3m+2)}
```
```math
d_0 = b_0 \qquad d_m = - \frac{d_{m-1}}{(3m-2)3m}
```

### Asymptotic Solution

For ``h₁`` on ``-2π/3 < \arg z < 4π/3``:
```math
h₁(z) ≈ α z^{-1/4} \exp(2/3 i z^{3/2} - 5πi/12) \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right)
```
```math
h₁'(z) ≈ α i z^{1/4} \exp(2/3 i z^{3/2} - 5πi/12) \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) \\
    - α/4 z^{-5/4} \exp(2/3 i z^{3/2} - 5πi/12) \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) \\
    - 3/2 α z^{-1/4} \exp(2/3 i z^{3/2} - 5πi/12) \left( \sum_{m=1} (-i)^m m C_m z^{-3m/2 - 1} \right)
```
where
```math
C_m = \frac{(9-4)(81-4)\cdots (9[2m-1]^2-4)}{2^{4m}3^m m!}
```
and
```math
h₁(z) ≈ h₁(z) + α z^{-1/4} \exp(-2/3 i z^{3/2} - 11πi/12) \left( 1 + \sum_{m_1} (i)^m C_m z^{-3m/2} \right)
```
```math
h₁'(z) ≈ h₁'(z) - α z^{1/4} \exp(-2/3 i z^{3/2} - 11πi/12) \left( 1 + \sum_{m_1} (i)^m C_m z^{-3m/2} \right) \\
    - α/4 z^{-5/4} \exp(-2/3 i z^{3/2} - 11πi/12) \left( 1 + \sum_{m_1} (i)^m C_m z^{-3m/2} \right) \\
    - 3/2 α z^{-1/4} \exp(-2/3 i z^{3/2} - 11πi/12) \left( \sum_{m=1} i^m m C_m z^{-3m/2 - 1} \right)
```
for ``-4π/3 < \arg z < 0``.

And for ``h₂`` on ``-4π/3 < \arg z < 2π/3``:
```math
h₂(z) ≈ α z^{-1/4} \exp(-2/3 i z^{3/2} + 5πi/12) \left( 1 + \sum_{m=1} (i)^m C_m z^{-3m/2} \right)
```
```math
h₂'(z) ≈ -α z^{1/4} \exp(-2/3 i z^{3/2} + 5πi/12) \left( 1 + \sum_{m=1} (i)^m C_m z^{-3m/2} \right) \\
    - α/4 z^{-5/4} \exp(-2/3 i z^{3/2} + 5πi/12) \left( 1 + \sum_{m=1} (i)^m C_m z^{-3m/2} \right) \\
    - 3/2 α z^{-1/4} \exp(-2/3 i z^{3/2} + 5πi/12) \left( i^m m C_m z^{-3m/2 - 1} \right)
```
and
```math
h₂(z) ≈ h₂(z) + α z^{-1/4} \exp(2/3 i z^{3/2} + 11πi/12) \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right)
```
```math
h₂'(z) ≈ h₂'(z) + α z^{1/4} \exp(2/3 i z^{3/2} + 11πi/12) \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) \\
    - α/4 z^{-5/4} \exp(2/3 i z^{3/2} + 11πi/12) \left( 1 + \sum_{m=1} (-i)^m C_m z^{-3m/2} \right) \\
    - 3/2 α z^{-1/4} \exp(2/3 i z^{3/2} + 11πi/12) \left( \sum_{m=1} (-i)^m m C_m z^{-3m/2 - 1} \right)
```
for ``0 < \arg z < 4π/3``.

The coefficient ``\alpha = 2^{1/3} 3^{1/6} \pi^{-1/2}``.

## References

[^SCL1945]:

    The Staff of the Computation Library (1945), *Tables of the modified Hankel
    function of order one-third and of their derivatives.* Cambridge, MA: Harvard
    University Press.
