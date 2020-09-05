# Performance Tips
There are a number of ways to improve the performance of the solver given a particular problem.
If you are not satisfied with the performance there are a number of things you have to determine first. Is the solver slow because
- it's the first time you ran it in the current Julia session or
- is it because the solver needs a lot of iterations (convergence) or
- is each iteration or the initial factorisation slow (computational performance)?

Let's see how each point can be addressed.

## First run
Whenever a new Julia session is started, the first run will trigger a compilation of all functions based on their arguments used in your script. This will slow the first execution of COSMO down. After that Julia will call the fast compiled functions. To get around this, you can either keep your current Julia session open and discard the first run. Alternatively, if your problem is very large, you could solve a small version of your problem first to trigger the compilation. Another option is to use [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl) to save compiled functions into a *sysimage* that can be loaded at startup.

## Solver Timings
It is often instructive to look at the detailed solver timing results for your problem. This can reveal where most of the time is spent. To achieve this, run COSMO with the setting `verbose_timing = true`. After solving the problem with `result = COSMO.optimize!(model)` you can look at `result.times` for a breakdown of the times spent in different parts of the algorithm, see [Timings](@ref) for more details. Especially take a look at the ratio of factorisation time and iteration time. If you use `JuMP` to solve the problem, you can take a look at the timings with `backend(model).optimizer.model.optimizer.results.times`.

## Convergence
It is possible that COSMO converges slowly, i.e. needs a large number of iterations, for your problem given its default parameters.

### Parameters
You could try changing any of the following parameters:
- `rho`: The initial algorithm step parameter has a large influence on the convergence. Try different values between `1e-5` and `10`.
- `adaptive_rho = false`: You can try to disable the automatic rho adaption and use different rho values.
- `adaptive_rho_interval`: This specifies after how many iterations COSMO tries to adapt the rho parameter. You can also set `adaptive_rho_interval = 0` which adapts the rho parameter after the time spent iterating passes 40% of the factorisation time. This is currently the default in [OSQP](https://osqp.org/) and works well with QPs.
- `alpha = 1.0`: This disables the over-relaxation that is used in the algorithm. We recommend values between `1.0 - 1.6`.
- `scaling = 0`: This disables the problem scaling.
- `eps_abs` and `eps_rel`: Check the impact of modifying the stopping accuracies.

### Use warm starting
The number of iterations can be dramatically decreased by providing a good initial guess for `x`, `s` and `y`. Examples where warm starting is commonly used are model predictive control and portfolio backtests, see [Warm starting](@ref).

## Computational performance
If the convergence of the algorithm is not an issue, there are still a number of steps you can take to make COSMO faster.

### Intel MKL BLAS/LAPACK
We experienced significant performance improvements on Intel CPUs if Julia is compiled with MKL BLAS. This is because Julia's linear algebra function will use Intel MKL BLAS and LAPACK functions that are optimised for Intel hardware. The effect is especially significant for SDPs because most of the time is spent in the LAPACK function `syevr`. If you are running Julia on Intel hardware, an easy way to compile Julia with MKL is to add and build the MKL package, see [MKL.jl](https://github.com/JuliaComputing/MKL.jl). To verify your current BLAS vendor you can use `julia> LinearAlgebra.BLAS.vendor()`.

### Linear system solver
COSMO uses `QDLDL.jl` as the default linear system solver. In our experience this seems to be a competitive choice until about `1e5 - 1e6` nonzeros in the constraint matrix. After that it is worth trying one of the indirect system solvers, such as CG or MINRES. Furthermore, we also recommend trying Pardiso (or MKLPardiso) for problems of that dimension. More details can be found here: [Linear System Solver](@ref).

### Custom cones
In some cases the computations can be speed-up if certain constraints in the problem allow the implementation of a fast projection function. We allow the user to define their own custom convex cone with a corresponding projection function.
The custom cone has to be defined as `struct CustomCone{T} <: COSMO.AbstractConvexSet{T}`. Furthermore, the user has to define a function that projects an input vector `x` onto the custom cone, i.e. `function COSMO.project!(x::AbstractVector{T}, C::CustomCone{T}) where {T} ... end`.

### Multithreading
COSMO allows the execution of the projection step for multiple constraints in parallel using Julia's multithreading features. This is currently not enabled in the tagged release because of stability issues in earlier Julia versions. To use multithreading checkout the branch `with_multi_threading`, which we keep in sync with the latest tagged release. This can be installed via Julia's package manager with  `pkg> add COSMO#with_multi_threading`. Afterwards, before starting Julia, set `export JULIA_NUM_THREADS=[NUMBER_LOGICAL_CORES_HERE]`. In Julia you can verify the number of threads with `julia> Threads.nthreads()`.

Notice that the extra overhead for multithreading can slow the solver down if the problem is small. However, we noticed significant performance improvements if the problem contained multiple positive semidefinite constraints or when one large constraint was decomposed. In that case it also helps to restrict the number of BLAS threads per Julia thread with `julia> BLAS.set_num_threads(1)` to prevent oversubscription of the available cores.

Multithreading can also be used in the factorisation step if the Pardiso or MKLPardiso solver are selected. This is only advisable for constraint matrices with more than `1e5` nonzeros.

### Chordal decomposition and Clique merging
When solving large structured and sparse SDPs significant performance improvements are achievable if the problem is passed to COSMO in the right way. This means the solver has to be able to infer the structure of the positive semidefinite variable from the constraint. See the section on [Chordal Decomposition](@ref) for more details. In some cases the primal SDP doesn't allow decomposition but the dual SDP does, consider the [Maximum Cut Problem](@ref) and the [Relaxed Two-Way Partitioning Problem](@ref) for examples.

If the problem is decomposable it is also worth experimenting with different clique merging strategies to see how they impact the performance. More details can be found here: [Clique merging](@ref).
