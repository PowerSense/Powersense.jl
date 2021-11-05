# Powersense.jl

<img src="https://powersense.github.io//assets/Powersense_logo.png" align="left" width="200" alt="PowerSense logo">

<a href="https://www.powersense.io/)"><img src="https://img.shields.io/badge/powered%20by-PowerSENSE-blue"/></a>
![Run tests](https://github.com/PowerSense/Powersense.jl/workflows/Run%20tests/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/PowerSense/Powersense.jl/branch/master/graph/badge.svg?token=SUH4VPE41D)](https://codecov.io/gh/PowerSense/Powersense.jl)
<!-- [![Documentation](https://github.com/PowerSense/Powersense.jl/workflows/Documentation/badge.svg)](https://www.powersense.io/) -->

The package has the following features implemented and ready to use.

- `OPT`: The nonlinear programming optimization component implements sequential linear programming method for continuous nonlinear optimization. The package currently implements a line search algorithm.
- `OPF`: The AC optimal power flow (AC-OPF) component implements 7 different AC-OPF formulations. These formulations have different sparsities and are a combination of varieties of different approaches of modeling voltages, admittance matrix, and branch flows. The formulations can be solved using the `OPT` feature of the `Powersense.jl` or using an external NLP solver.


## Installation

```julia
import Pkg
Pkg.add("Powersense")
```



## NLP Solver Example

Consider the following quadratic optimization problem

```julia
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

## OPF formulations Example

Consider solving PNPAPVmodel ACOPF formulations. This can be solved using the following code snippet:
```julia
# Load packages
using Powersense, Ipopt

# Build Powersense Data model. Path is the address to where PSSE or MATPOWER file types
Data = create_PowersenseData(path)

run_opf!(Data, solver = Ipopt.Optimizer, obj_type = "linear", formulation = PNPAPVmodel);
```

# Acknowledgements

The package is part of the [PowerSense Lab](https://www.powersense.io/) which is owned and maintained by [Sayed Abdullah Sadat](https://www.sayedsadat.com/).
