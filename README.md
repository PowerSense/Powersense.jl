# Powersense.jl
![Run tests](https://github.com/exanauts/ActiveSetMethods/workflows/Run%20tests/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/exanauts/ActiveSetMethods/branch/master/graph/badge.svg)](https://codecov.io/gh/exanauts/ActiveSetMethods)

This is a Julia package that implements active set methods for continuous nonlinear optimization.
The package currently implements a sequential linear programming method based on line search.

## Installation

```julia
]add https://github.com/exanauts/ActiveSetMethods.jl
```

## Example

Consider the following quadratic optimization problem

```
min   x^2 + x 
s.t.  x^2 - x = 2
```

This problem can be solved by the following code snippet:
```julia
# Load packages
using ActiveSetMethods, JuMP
using GLPK # can be any LP solver

# Number of variables
n = 1

# Build nonlinear problem model via JuMP
model = Model(optimizer_with_attributes(ActiveSetMethods.Optimizer, "external_optimizer" => GLPK.Optimizer))
@variable(model, x)
@objective(model, Min, x^2 + x)
@NLconstraint(model, x^2 - x == 2)

# Solve optimization problem with Nlopt
JuMP.optimize!(model)

# Retrieve solution
Xsol = JuMP.value.(X)
```

## Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.
