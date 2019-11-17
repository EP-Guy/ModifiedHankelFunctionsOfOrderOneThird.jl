# Modified Hankel Functions Of Order One Third and Their Derivatives

[![Build Status](https://travis-ci.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl.svg?branch=master)](https://travis-ci.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl) [![Build status](https://ci.appveyor.com/api/projects/status/w115vkl46t4nj4ui?svg=true)](https://ci.appveyor.com/project/EP-Guy/modifiedhankelfunctionsoforderonethird) [![DOI](https://zenodo.org/badge/156012814.svg)](https://zenodo.org/badge/latestdoi/156012814)


Solutions to Stokes' differential equation:

![\frac{\mathrm{d}^2u}{\mathrm{d}z^2} + zu = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cmathrm%7Bd%7D%5E2u%7D%7B%5Cmathrm%7Bd%7Dz%5E2%7D%20&plus;%20zu%20%3D%200)

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

![h_1(z) = \frac{k}{i\pi} \int_{L_1} e^{zt + \frac{t^3}{3}} \,\mathrm{d}t](https://latex.codecogs.com/svg.latex?h_1%28z%29%20%3D%20%5Cfrac%7Bk%7D%7Bi%5Cpi%7D%20%5Cint_%7BL_1%7D%20e%5E%7Bzt%20&plus;%20%5Cfrac%7Bt%5E3%7D%7B3%7D%7D%20%5C%2C%5Cmathrm%7Bd%7Dt)

and

![h_2(z) = \frac{k^*}{-i\pi} \int_{L_2} e^{zt + \frac{t^3}{3}} \,\mathrm{d}t](https://latex.codecogs.com/svg.latex?h_2%28z%29%20%3D%20%5Cfrac%7Bk%5E*%7D%7B-i%5Cpi%7D%20%5Cint_%7BL_2%7D%20e%5E%7Bzt%20&plus;%20%5Cfrac%7Bt%5E3%7D%7B3%7D%7D%20%5C%2C%5Cmathrm%7Bd%7Dt)

where ![k = (12)^\frac{1}{6} e^{\left(-\frac{\pi i}{6} \right)}](https://latex.codecogs.com/svg.latex?k%20%3D%20%2812%29%5E%5Cfrac%7B1%7D%7B6%7D%20e%5E%7B%5Cleft%28-%5Cfrac%7B%5Cpi%20i%7D%7B6%7D%20%5Cright%29%7D)

The contours of integration _L₁_ and _L₂_ are

![contoursofintegration](contoursofintegration.svg)

with π/6 < _w_ < π/2. We take _w_ = π/3.

## Solutions

Two solution approaches are used. If `abs2(z) < 36`, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.

### Power series

Stokes' equation may be solved in a power series of _z_, valid in the entire complex plane,

![h_1(z) = g + \frac{i\sqrt{3}}{3} ( g - 2f), \quad h_2(z) = g - \frac{i\sqrt{3}}{3} ( g - 2f)](https://latex.codecogs.com/svg.latex?h_1%28z%29%20%3D%20g%20&plus;%20%5Cfrac%7Bi%5Csqrt%7B3%7D%7D%7B3%7D%20%28%20g%20-%202f%29%2C%20%5Cquad%20h_2%28z%29%20%3D%20g%20-%20%5Cfrac%7Bi%5Csqrt%7B3%7D%7D%7B3%7D%20%28%20g%20-%202f%29)

where

![f(z) = \frac{2^{1/3}}{\Gamma(\frac{2}{3})} \left[ 1 + \sum_1^\infty \frac{(-)^m (3m-2)(3m-5)\cdots 4\cdot 1}{(3m)!} z^{3m} \right]](https://latex.codecogs.com/svg.latex?f%28z%29%20%3D%20%5Cfrac%7B2%5E%7B1/3%7D%7D%7B%5CGamma%28%5Cfrac%7B2%7D%7B3%7D%29%7D%20%5Cleft%5B%201%20&plus;%20%5Csum_1%5E%5Cinfty%20%5Cfrac%7B%28-%29%5Em%20%283m-2%29%283m-5%29%5Ccdots%204%5Ccdot%201%7D%7B%283m%29%21%7D%20z%5E%7B3m%7D%20%5Cright%5D)

![g(z) = \frac{2^{1/3}}{3^{2/3}\Gamma(\frac{4}{3})} \left[ z + \sum_1^\infty \frac{(-)^m (3m-1)(3m-4)\cdots 5\cdot 2}{(3m+1)!} z^{3m+1} \right]](https://latex.codecogs.com/svg.latex?g%28z%29%20%3D%20%5Cfrac%7B2%5E%7B1/3%7D%7D%7B3%5E%7B2/3%7D%5CGamma%28%5Cfrac%7B4%7D%7B3%7D%29%7D%20%5Cleft%5B%20z%20&plus;%20%5Csum_1%5E%5Cinfty%20%5Cfrac%7B%28-%29%5Em%20%283m-1%29%283m-4%29%5Ccdots%205%5Ccdot%202%7D%7B%283m&plus;1%29%21%7D%20z%5E%7B3m&plus;1%7D%20%5Cright%5D)

### Asymptotic expansion

The asymptotic expansions can be used to estimate _h₁_, _h₂_, and their derivatives, although in general with less accuracy than the power series. Two expansions are required depending on the value of `arg z`. The existence of two expressions of different forms which represent asymptotically the same integral function is an example of [Stokes' phenomenon](https://en.wikipedia.org/wiki/Stokes_phenomenon).

The expansion for _h₁_ for -2π/3 < arg _z_ < 4π/3 is

![h_1(z) \approx \alpha z^{-1/4} e^{\frac{2}{3}iz^{3/2} - \frac{5\pi i}{12}} \left[ 1 + \sum_{m=1} (-i)^m C_m z^\frac{-3m}{2} \right]](https://latex.codecogs.com/svg.latex?h_1%28z%29%20%5Capprox%20%5Calpha%20z%5E%7B-1/4%7D%20e%5E%7B%5Cfrac%7B2%7D%7B3%7Diz%5E%7B3/2%7D%20-%20%5Cfrac%7B5%5Cpi%20i%7D%7B12%7D%7D%20%5Cleft%5B%201%20&plus;%20%5Csum_%7Bm%3D1%7D%20%28-i%29%5Em%20C_m%20z%5E%5Cfrac%7B-3m%7D%7B2%7D%20%5Cright%5D)

where

![C_m = \frac{(9-4)(81-4)\cdots (9[2m-1]^2-4)}{2^{4m}3^m m!}](https://latex.codecogs.com/svg.latex?C_m%20%3D%20%5Cfrac%7B%289-4%29%2881-4%29%5Ccdots%20%289%5B2m-1%5D%5E2-4%29%7D%7B2%5E%7B4m%7D3%5Em%20m%21%7D)

See the source for the full sets of solutions.

## References

The Staff of the Computation Library (1945), _Tables of the modified Hankel function of order one-third and of their derivatives_. Cambridge, MA: Harvard University Press.

## Citing

We encourage you to cite this package if used in scientific work. See the Zenodo
badge above or refer to [CITATION.bib](CITATION.bib).
