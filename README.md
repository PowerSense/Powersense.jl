# Powersense.jl

<img src="https://powersense.github.io//assets/Powersense_logo.png" align="left" width="200" alt="PowerSense logo">

[![codecov](https://img.shields.io/badge/powered%20by-PowerSENSE-blue)](https://www.powersense.io/)
![Run tests](https://github.com/PowerSense/Powersense.jl/workflows/Run%20tests/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/PowerSense/Powersense.jl/branch/master/graph/badge.svg?token=SUH4VPE41D)](https://codecov.io/gh/PowerSense/Powersense.jl)
[![Documentation](https://github.com/PowerSense/Powersense.jl/workflows/Documentation/badge.svg)](https://www.powersense.io/)

The NLP solver of this Julia package implements active set methods for continuous nonlinear optimization. The package currently implements a sequential linear programming method based on line search.


## Installation

```julia
import Pkg
Pkg.add("Powersense")
```



## NLP Solver Example

Consider the following quadratic optimization problem

```
min   x^2 + x 
s.t.  x^2 - x = 2
```

This problem can be solved by the following code snippet:
```julia
# Load packages
using Powersense, JuMP
using GLPK # can be any LP solver

# Number of variables
n = 1

# Build nonlinear problem model via JuMP
model = Model(optimizer_with_attributes(Powersense.Optimizer, "external_optimizer" => GLPK.Optimizer))
@variable(model, x)
@objective(model, Min, x^2 + x)
@NLconstraint(model, x^2 - x == 2)

# Solve optimization problem with Nlopt
JuMP.optimize!(model)

# Retrieve solution
Xsol = JuMP.value.(X)
```

# Acknowledgements

The package is part of the [PowerSense Lab](https://www.powersense.io/) which is owned and maintained by [Sayed Abdullah Sadat](https://www.sayedsadat.com/).
