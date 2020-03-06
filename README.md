# Modified Hankel Functions Of Order One Third and Their Derivatives

[![Build Status](https://travis-ci.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl.svg?branch=master)](https://travis-ci.com/fgasdia/ModifiedHankelFunctionsOfOrderOneThird.jl) [![Build status](https://ci.appveyor.com/api/projects/status/w115vkl46t4nj4ui?svg=true)](https://ci.appveyor.com/project/EP-Guy/modifiedhankelfunctionsoforderonethird) [![DOI](https://zenodo.org/badge/156012814.svg)](https://zenodo.org/badge/latestdoi/156012814) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://fgasdia.github.io/ModifiedHankelFunctionsOfOrderOneThird.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://fgasdia.github.io/ModifiedHankelFunctionsOfOrderOneThird.jl/dev)


Solutions to Stokes' differential equation:

<p align="center"><img src="/tex/7a9703279d8af8b6c816345789cfc3d6.svg?invert_in_darkmode&sanitize=true" align=middle width=95.89569494999999pt height=35.77743345pt/></p>

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

An independent pair of solutions, valid for all values of <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/>, is

<p align="center"><img src="/tex/0881e9f68f46ac526be6ce0622c73547.svg?invert_in_darkmode&sanitize=true" align=middle width=173.4397896pt height=39.1573875pt/></p>

and

<p align="center"><img src="/tex/43f1a610ddc4c7f2989ac468c8038443.svg?invert_in_darkmode&sanitize=true" align=middle width=186.22521554999997pt height=39.1573875pt/></p>

where <img src="/tex/6a1370aa3b609ee19cdf955664a9c204.svg?invert_in_darkmode&sanitize=true" align=middle width=116.74815569999998pt height=36.4155132pt/>.

The contours of integration <img src="/tex/929ed909014029a206f344a28aa47d15.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73978854999999pt height=22.465723500000017pt/> and <img src="/tex/4327ea69d9c5edcc8ddaf24f1d5b47e4.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73978854999999pt height=22.465723500000017pt/> are

![contoursofintegration](contoursofintegration.svg)

with <img src="/tex/8f8e464d2606edffa767ab489706c317.svg?invert_in_darkmode&sanitize=true" align=middle width=108.84312119999998pt height=24.65753399999998pt/>. We take <img src="/tex/7a89758c6f0f5b13c3254cb7dfd3b713.svg?invert_in_darkmode&sanitize=true" align=middle width=60.52698299999998pt height=24.65753399999998pt/>.

## Solutions

Two solution approaches are used. If `abs2(z) < 36`, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.

### Power series

Stokes' equation may be solved in a power series of <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/>, valid in the entire complex plane,

<p align="center"><img src="/tex/2386478755dc5db6621d0ba4086a70b8.svg?invert_in_darkmode&sanitize=true" align=middle width=382.35967934999996pt height=36.65224035pt/></p>

where

<p align="center"><img src="/tex/c957005964fab87d65f9dfb65fd40d2e.svg?invert_in_darkmode&sanitize=true" align=middle width=416.97635099999997pt height=49.315569599999996pt/></p>

<p align="center"><img src="/tex/7538b8a402fd0a8b8104ce9d3a82a397.svg?invert_in_darkmode&sanitize=true" align=middle width=461.26299285pt height=49.315569599999996pt/></p>

### Asymptotic expansion

The asymptotic expansions can be used to estimate <img src="/tex/5a95dbebd5e79e850a576db54f501ab8.svg?invert_in_darkmode&sanitize=true" align=middle width=16.02366149999999pt height=22.831056599999986pt/>, <img src="/tex/0f7cea0b89929faf20eda59174bc247f.svg?invert_in_darkmode&sanitize=true" align=middle width=16.02366149999999pt height=22.831056599999986pt/>, and their derivatives, although in general with less accuracy than the power series. Two expansions are required depending on the value of `arg z`. The existence of two expressions of different forms which represent asymptotically the same integral function is an example of [Stokes' phenomenon](https://en.wikipedia.org/wiki/Stokes_phenomenon).

The expansion for <img src="/tex/5a95dbebd5e79e850a576db54f501ab8.svg?invert_in_darkmode&sanitize=true" align=middle width=16.02366149999999pt height=22.831056599999986pt/> for <img src="/tex/2bc9dc1ce765d768fd71a26518a924af.svg?invert_in_darkmode&sanitize=true" align=middle width=160.06845525pt height=24.65753399999998pt/> is

<p align="center"><img src="/tex/e7e4ce30240a2c855f796c7ae9e72946.svg?invert_in_darkmode&sanitize=true" align=middle width=372.81839924999997pt height=49.315569599999996pt/></p>

where

<p align="center"><img src="/tex/c23b4cbc617449d1f0f90b637fc3a962.svg?invert_in_darkmode&sanitize=true" align=middle width=296.41710945pt height=35.77743345pt/></p>

See the source for the full sets of solutions.

## References

The Staff of the Computation Library (1945), _Tables of the modified Hankel function of order one-third and of their derivatives_. Cambridge, MA: Harvard University Press.

## Citing

We encourage you to cite this package if used in scientific work. See the Zenodo
badge above or refer to [CITATION.bib](CITATION.bib).
