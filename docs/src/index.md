# Modified Hankel Functions of Order One Third

Solutions to Stokes' differential equation:

![\frac{\mathrm{d}^2u}{\mathrm{d}z^2} + zu = 0](https://latex.codecogs.com/svg.latex?%5Cfrac%7B%5Cmathrm%7Bd%7D%5E2u%7D%7B%5Cmathrm%7Bd%7Dz%5E2%7D%20&plus;%20zu%20%3D%200)

```@contents
```

## Solutions

Two solution approaches are used. If `abs2(z) < 36`, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.

## Usage

```julia
using ModifiedHankelFunctionsOfOrderOneThird

h1, h2, h1′, h2′ = modifiedhankel(z)
```
