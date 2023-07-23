# Linear System Solver

One major step of COSMO's ADMM algorithm is solving a linear system of equations at each iteration. Fortunately, the left-hand matrix is only dependent on the problem data and therefore only needs to be factored once. Depending on the problem class this factorisation can be the computationally most expensive step of the algorithm (LPs, QPs). See the [Method](@ref) section for a more detailed description of the linear system.

COSMO allows you to specify the linear system solver that performs the factorisation and back-substitution. We also support indirect system solver which are useful for very large problems where a factorisation becomes inefficient. The table below shows the currently supported linear system solver:



Type | Solver | Description
--- | --- | ---
CholmodKKTSolver | Cholmod | Julia's default linear system solver (from SuiteSparse)
QdldlKKTSolver | QDLDL | For more information [QDLDL.jl](https://github.com/oxfordcontrol/QDLDL.jl)
PardisoDirectKKTSolver | Pardiso (direct) | Pardiso 6.0 direct solver
PardisoIndirectKKTSolver | Pardiso (indirect) | Pardiso 6.0 indirect solver
MKLPardisoKKTSolver | Intel MKL Pardiso | Pardiso optimised for Intel platforms
CGIndirectKKTSolver | IterativeSolvers.jl | Conjugate Gradients on the reduced KKT linear system.
MINRESIndirectKKTSolver | IterativeSolvers.jl | MINRES on the (full) KKT linear system.

!!! note
    To use the Pardiso and Intel MKL Pardiso solver, you have to install the respective libraries and the corresponding Julia wrapper. For more information about installing these, visit the [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl) repository page. Likewise in order to use Indirect(Reduced)KKTSolver you have to install [IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl) (v0.9+) and [LinearMaps.jl](https://github.com/Jutho/LinearMaps.jl). We are using the `Requires` package for lazy loading of code related to `Pardiso` and `IterativeSolvers`. This means in order to use `Pardiso` / `IterativeSolvers`, you'll have to load these packages alongside `COSMO`, i.e. `using Pardiso` and `using IterativeSolvers, LinearMaps`.

COSMO uses the `Cholmod` linear system solver by default. You can specify a different solver in the settings by using the `kkt_solver` keyword and the respective type:

```julia
settings = COSMO.Settings(kkt_solver = CholmodKKTSolver)

```

COSMO also allows you to pass in solver-specific options with the `with_options(solver_type, args...; kwargs...)` syntax. For example, if you want to use Pardiso with verbose printing use the following syntax:
```julia
settings = COSMO.Settings(kkt_solver = with_options(PardisoDirectKKTSolver, msg_level_on = true))
```

Likewise, `CGIndirectKKTSolver` and `MINRESIndirectKKTSolver` are also parameterizable with `with_options(solver_type, args...; kwargs...)` and accept the following arguments:

 Argument | Description | Values (default)
-------------- |   :-------------- |   :--------------
`tol_constant::T` and `tol_exponent::T` | Parameter that defines the solution tolerance for the iterative solvers across iterations. In particular, the solution tolerance at every iteration is defined as ``\text{tol\_constant} /  \text{iteration}^{\text{tol\_exponent}}`` | 1.0, 1.5


This also works if you want to use this configuration with JuMP:

```julia
model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "kkt_solver" => with_options(PardisoDirectKKTSolver, msg_level_on = true));

```
Or alternatively:
```julia
model = JuMP.Model(COSMO.Optimizer);
set_optimizer_attribute(model, "kkt_solver", with_options(PardisoDirectKKTSolver, msg_level_on = true));
```