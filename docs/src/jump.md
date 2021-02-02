# JuMP Interface
Our [JuMP](https://github.com/JuliaOpt/JuMP.jl/) interface allows you to describe and modify your optimisation problem with JuMP and use COSMO as the backend solver. The interface is defined in `/src/MOIWrapper.jl`.

!!! note
    We assume here that the latest JuMP release (`~0.21.0`) is used.

## Use COSMO
To specify COSMO as the solver for your JuMP model, load the solver module with `using COSMO` and then pass a `COSMO.Optimizer` when initialising the JuMP model:
```julia
m = JuMP.Model(COSMO.Optimizer);
```

## Specify Solver Settings
Solver-specific settings can be passed using the `optimizer_with_attributes()` function. For example, if you want to adjust the maximum number of iterations and turn on verbose printing use
```julia
m = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "max_iter" => 5000, "verbose" => true));
```
Note that the attributes are passed as key-value pairs and the keys are strings. This is slightly different to using the native COSMO interface. Equivalently, one can also use:
```julia
m = JuMP.Model(COSMO.Optimizer);
set_optimizer_attribute(m, "max_iter", 5000)
set_optimizer_attribute(m, "verbose", true)
```

The full list of available settings can be found in the [Settings](#settings) section. All of them should be compatible with JuMP.

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
To query the number of iterations of the algorithm use:
```julia
iter = MOI.get(m, COSMO.ADMMIterations())
```


For more information on how to use JuMP check the [JuMP documentation](http://www.juliaopt.org/JuMP.jl/stable/).
