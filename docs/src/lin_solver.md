# Linear System Solver

One major step of COSMO's ADMM algorithm is solving a linear system of equations at each iteration. Fortunately, the left-hand matrix is only dependent on the problem data and therefore only needs to be factored once. Depending on the problem class this factorisation can be the computationally most expensive step of the algorithm (LPs, QPs). See the [Method](@ref) section for a more detailed description of the linear system.

COSMO allows you to specify the linear system solver that performs the factorisation and back-substitution. We also support indirect system solver which are useful for very large problems where a factorisation becomes inefficient. The table below shows the currently supported linear system solver:



Type | Solver | Description
--- | --- | ---
QDLDLKKTSolver | QDLDL | For more information [QDLDL.jl](https://github.com/oxfordcontrol/QDLDL.jl)
CholmodKKTSolver | Cholmod | Julia's default linear system solver (by SuiteSparse)
PardisoDirectKKTSolver | Pardiso (direct) | Pardiso 6.0 direct solver
PardisoIndirectKKTSolver | Pardiso (indirect) | Pardiso 6.0 indirect solver
MKLPardisoKKTSolver | Intel MKL Pardiso | Pardiso optimised for Intel platforms

!!! note
    To use the Pardiso and Intel MKL Pardiso solver, you have to install the respective libraries and the corresponding Julia wrapper. For more information about installing these, visit the [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl) repository page.

COSMO uses the QDLDL linear system solver by default. You can specify a different solver in the settings by using the `kkt_solver` keyword and the respective type:

```julia
settings = COSMO.Settings(kkt_solver = CholmodKKTSolver)

```

COSMO also allows you to pass in solver-specific options. If you want to use Pardiso with verbose printing use the `with_options(solver_type, args...; kwargs...)` syntax:


```julia
settings = COSMO.Settings(kkt_solver = with_options(PardisoDirectKKTSolver, msg_level_on = true))

```

This also works if you want to use this configuration with JuMP:


```julia
model = JuMP.Model(with_optimizer(COSMO.Optimizer, kkt_solver = with_options(PardisoDirectKKTSolver, msg_level_on = true));

```