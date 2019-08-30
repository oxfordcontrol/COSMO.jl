# JuMP Interface
Our [JuMP](https://github.com/JuliaOpt/JuMP.jl/) interface allows you to describe and modify your optimisation problem with JuMP and use COSMO as the backend solver. The interface is defined in `/src/MOIWrapper.jl`.

!!! note
    COSMO requires the newest JuMP release that is based on the MathOptInterface package.

## Use COSMO
To specify COSMO as the solver for your JuMP model, load the solver module with `using COSMO` and then use the `with_optimizer()` function when initialising the JuMP model:
```julia
m = JuMP.Model(with_optimizer(COSMO.Optimizer);
```

## Specify Solver Settings
Solver-specific settings can be passed after the `COSMO.Optimizer` object. For example, if you want to adjust the maximum number of iterations and turn on verbose printing use
```julia
m = JuMP.Model(with_optimizer(COSMO.Optimizer, max_iter = 5000, verbose = true);
```
The full list of available settings can be found in the [Settings](#settings) section.

## Results
After solving the problem the result can be obtained using the standard JuMP commands. To see if the optimisation was successful use
```julia
JuMP.termination_status(m)
JuMP.primal_status(m)
```
If a solution is available, the optimal objective value can be queried using
```julia
JuMP.objective_value(m)
```
and the value of a decision variable `x` can be obtained with
```julia
JuMP.value.(x)
```
For more information on how to use JuMP check the [JuMP documentation](http://www.juliaopt.org/JuMP.jl/stable/).