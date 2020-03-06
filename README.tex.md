# Modified Hankel Functions Of Order One Third and Their Derivatives

[![Build Status](https://travis-ci.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl.svg?branch=master)](https://travis-ci.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl) [![Build status](https://ci.appveyor.com/api/projects/status/w115vkl46t4nj4ui?svg=true)](https://ci.appveyor.com/project/EP-Guy/modifiedhankelfunctionsoforderonethird) [![DOI](https://zenodo.org/badge/156012814.svg)](https://zenodo.org/badge/latestdoi/156012814)


Solutions to Stokes' differential equation:

$$ \frac{\mathrm{d}^2u}{\mathrm{d}z^2} + zu = 0 $$

From _Tables of the Modified Hankel Functions of Order One-Third and of their Derivatives_:

> Its only singularity is an irregular singularity at infinity. The equation occurs in the description of simple cases of diffraction and of refraction of waves.
> The general solution of [Stokes' equation] can be written in terms of Bessel functions of order one-third. The tabulation of these Bessel functions for complex arguments would make possible the computation of solutions of [Stokes' equation] for complex arguments. The direct tabulation of solutions of [Stokes' equation] should, however, be preferred to that of Bessel functions of order one-third. Unlike Bessel's equation, **[Stokes' equation] has no singularity in the finite complex plane and its solutions are single-valued**, whereas the Bessel functions of order one-third are not.

## Usage

To install:
```julia
]add https://github.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl
```

Then:
```julia
using ModifiedHankelFunctionsOfOrderOneThird

h1, h2, h1′, h2′ = modifiedhankel(z)
```

## The functions _h₁_ and _h₂_

An independent pair of solutions, valid for all values of _z_, is

$$ h_1(z) = \frac{k}{i\pi} \int_{L_1} e^{zt + \frac{t^3}{3}} \,\mathrm{d}t $$

and

$$ h_2(z) = \frac{k^\ast}{-i\pi} \int_{L_2} e^{zt + \frac{t^3}{3}} \,\mathrm{d}t $$

where $k = (12)^\frac{1}{6} e^{\left(-\frac{\pi i}{6} \right)}$

The contours of integration _L₁_ and _L₂_ are

![contoursofintegration](contoursofintegration.svg)

with π/6 < _w_ < π/2. We take _w_ = π/3.

## Solutions

Two solution approaches are used. If `abs2(z) < 36`, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.

### Power series

Stokes' equation may be solved in a power series of _z_, valid in the entire complex plane,

$$ h_1(z) = g + \frac{i\sqrt{3}}{3} ( g - 2f), \quad h_2(z) = g - \frac{i\sqrt{3}}{3} ( g - 2f) $$

where

$$ f(z) = \frac{2^{1/3}}{\Gamma(\frac{2}{3})} \left[ 1 + \sum_1^\infty \frac{(-)^m (3m-2)(3m-5)\cdots 4\cdot 1}{(3m)!} z^{3m} \right] $$

$$ g(z) = \frac{2^{1/3}}{3^{2/3}\Gamma(\frac{4}{3})} \left[ z + \sum_1^\infty \frac{(-)^m (3m-1)(3m-4)\cdots 5\cdot 2}{(3m+1)!} z^{3m+1} \right] $$

### Asymptotic expansion

The asymptotic expansions can be used to estimate _h₁_, _h₂_, and their derivatives, although in general with less accuracy than the power series. Two expansions are required depending on the value of `arg z`. The existence of two expressions of different forms which represent asymptotically the same integral function is an example of [Stokes' phenomenon](https://en.wikipedia.org/wiki/Stokes_phenomenon).

The expansion for _h₁_ for -2π/3 < arg _z_ < 4π/3 is

$$ h_1(z) \approx \alpha z^{-1/4} e^{\frac{2}{3}iz^{3/2} - \frac{5\pi i}{12}} \left[ 1 + \sum_{m=1} (-i)^m C_m z^\frac{-3m}{2} \right] $$

where

$$ C_m = \frac{(9-4)(81-4)\cdots (9[2m-1]^2-4)}{2^{4m}3^m m!} $$

See the source for the full sets of solutions.

## References

The Staff of the Computation Library (1945), _Tables of the modified Hankel function of order one-third and of their derivatives_. Cambridge, MA: Harvard University Press.

## Citing

We encourage you to cite this package if used in scientific work. See the Zenodo
badge above or refer to [CITATION.bib](CITATION.bib).
