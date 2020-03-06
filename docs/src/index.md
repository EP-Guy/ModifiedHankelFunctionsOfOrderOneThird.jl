# Modified Hankel Functions of Order One Third

_Calculate modified Hankel functions of order one-third and their derivatives._

These special functions are solutions to Stokes' differential equation:
``\frac{\mathrm{d}^2u}{\mathrm{d}z^2} + zu = 0``

The modified Hankel functions are linearly related to Airy functions but are not
as well known. One application of the modified Hankel functions of order one-third
is as a solution to the wave equation for electromagnetic fields between
the boundaries of the [earth-ionosphere waveguide](https://en.wikipedia.org/wiki/Earth%E2%80%93ionosphere_waveguide).

## Solutions

Two solution approaches are used, as presented in [^SCL1945]. If `abs2(z) < 36`,
a power series solution is used. Otherwise, an asymptotic expansion is performed
because of floating point limits in the power series.

## Usage

```julia
using ModifiedHankelFunctionsOfOrderOneThird

h1, h2, h1′, h2′ = modifiedhankel(z)
```

[^SCL1945]:
  The Staff of the Computation Library (1945), *Tables of the modified Hankel function of
  order one-third and of their derivatives.* Cambridge, MA: Harvard University Press.
