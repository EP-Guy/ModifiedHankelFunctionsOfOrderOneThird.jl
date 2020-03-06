# Modified Hankel Functions of Order One Third

Solutions to Stokes' differential equation:

<p align="center"><img src="/docs/src/tex/7a9703279d8af8b6c816345789cfc3d6.svg?invert_in_darkmode&sanitize=true" align=middle width=95.89569494999999pt height=35.77743345pt/></p>

```@contents
```

## Solutions

Two solution approaches are used. If `abs2(z) < 36`, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.

## Usage

```julia
using ModifiedHankelFunctionsOfOrderOneThird

h1, h2, h1′, h2′ = modifiedhankel(z)
```
