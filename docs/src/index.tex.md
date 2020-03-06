# Modified Hankel Functions of Order One Third

Solutions to Stokes' differential equation:

$$ \frac{\mathrm{d}^2u}{\mathrm{d}z^2} + zu = 0 $$

```@contents
```

## Solutions

Two solution approaches are used. If `abs2(z) < 36`, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.

## Usage

```julia
using ModifiedHankelFunctionsOfOrderOneThird

h1, h2, h1′, h2′ = modifiedhankel(z)
```
