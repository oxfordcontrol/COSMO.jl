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

## Feasibility of solution
In case `JuMP.termination_status(m)` is `MOI.ITERATION_LIMIT`,
`JuMP.primal_status(m)` (resp. `JuMP.dual_status(m)`) should be checked to determine
the feasibility status of the primal (resp. dual) solutions.
A primal (resp. dual) solution is considered feasible (resported as `MOI.FEASIBLE_POINT`) if
```math
r < \epsilon_{\text{abs}} + s\epsilon_{\text{rel}}
```
where
* ``r`` is the primal (resp. dual) residual `r_prim` ``= \lVert Ax + s - b \rVert_\infty``
(resp. `r_dual` ``= \lVert P x + q - A^\top * μ \rVert_\infty``),
* ``s`` is the primal (resp. dual) scaling of this residual
`max_norm_prim` ``= \max(\lVert Ax \rVert_\infty, \lVert s \rVert_\infty, \lVert b \rVert_\infty)``
(resp. `max_norm_dual` ``= \max(\lVert P x \rVert_\infty, \lVert q \rVert_\infty, \lVert A^\top * μ \rVert_\infty)``).
* ``\epsilon_{\text{abs}}`` is the optimizer attribute `eps_abs` and
* ``\epsilon_{\text{rel}}`` is the optimizer attribute `eps_rel`.

A primal (resp. dual solution) is considered *nearly* feasible (resported as `MOI.NEARLY_FEASIBLE_POINT`)
if it is not feasible but the same inequality is satisfied where the right-hand side is multiplied by
the optimizer attribute `nearly_ratio`.

Otherwise, the solution is considered infeasible (resported as `MOI.INFEASIBLE_POINT`).
Note that this does not mean that the problem is infeasible, it only means that the current primal (resp. dual)
solution is infeasible.

Note that, in case `JuMP.termination_status(m)` is `MOI.ITERATION_LIMIT`, the feasibility of the solution
is checked with the current optimizer attributes, `eps_abs`, `eps_rel` and `nearly_ratio`.
This means that if they are changed after `optimize!`, the status may change!
This can be used to check whether the solution is feasible up to a relaxed tolerance.
For instance, the following could happen:
```julia
optimize!(m)
primal_status(m) # MOI.INFEASIBLE_POINT
# This means that `r_prim >= nearly_ratio * (eps_abs * max_norm_prim * eps_rel)`
set_optimizer_attribute(m, "eps_abs", 1e-4)
primal_status(m) # MOI.NEARLY_FEASIBLE_POINT
# This means that `1 <= r_prim / (eps_abs * max_norm_prim * eps_rel) < nearly_ratio`
```

The values of `r_prim`, `max_norm_prim`, `r_dual` and `max_norm_dual` can also be accessed
as the fields of the `res_info` struct that can be obtained as follows
```julia
res_info = MOI.get(m, COSMO.RawResult()).info
```
Then, the feasibility can either be checked manually as in
```julia
res_info.r_prim < eps_abs + res_info.max_norm_prim * eps_rel
```
or using one of the following:
```julia
COSMO.is_primal_feasible(res_info, eps_abs, eps_rel)
COSMO.is_primal_nearly_feasible(res_info, eps_abs, eps_rel)
COSMO.is_dual_feasible(res_info, eps_abs, eps_rel)
COSMO.is_dual_nearly_feasible(res_info, eps_abs, eps_rel)
```


For more information on how to use JuMP check the [JuMP documentation](http://www.juliaopt.org/JuMP.jl/stable/).
