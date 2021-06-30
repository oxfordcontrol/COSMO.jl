# Acceleration 

By default COSMO's ADMM algorithm is wrapped in a safeguarded acceleration method to achieve faster convergence to higher precision. COSMO uses accelerators from the [COSMOAccelerators.jl](https://github.com/oxfordcontrol/COSMOAccelerators.jl) package.

By default, the solver uses the accelerator type `AndersonAccelerator{T, Type2{QRDecomp}, RestartedMemory, NoRegularizer}`. This is the classic type-II Anderson acceleration method where the least squares subproblem is solved using an updated QR method. Moreover, the method is restarted, i.e. the history of iterates is deleted, after `mem` steps and no regularisation for the least-squares method is used.

In addition, the method is safeguarded (`safeguard = true`), i.e. the residual-norm of the accelerated point can not deviate too much from the current point. Otherwise, the point is discarded and the ADMM algorithm performs a normal step instead.

The acceleration method can be altered as usual via the solver settings and the `accelerator` keyword. To deactivate acceleration pass an `EmptyAccelerator`:
```julia
settings = COSMO.Settings(accelerator = EmptyAccelerator)
```
To use the default accelerator but with a different memory size (number of stored iterates) use:
```julia
settings = COSMO.Settings(accelerator = with_options(AndersonAccelerator, mem = 15))
```
To turn the safeguarding off use:
```julia
settings = COSMO.Settings(safeguard = false)
```
To use an Anderson Accelerator of Type-I with a rolling-memory (oldest iterate replaced by newest) approach, use:
```julia
settings = COSMO.Settings(accelerator = AndersonAccelerator{Float64, Type1, RollingMemory, NoRegularizer})
```
For more fine-grained control look at the implementation of the accelerator [here](https://github.com/oxfordcontrol/COSMOAccelerators.jl/blob/master/src/anderson_accelerator.jl#L71).

When JuMP is used, the accelerator settings can be passed in the usual way:
```julia
model = JuMP.Model(optimizer_with_attributes(COSMO.Optimizer, "accelerator" => with_options(AndersonAccelerator, mem = 15)))
```
Or using the `set_optimizer_attribute()` method:
```julia
model = JuMP.Model(COSMO.Optimizer);
set_optimizer_attribute(model, "accelerator", with_options(AndersonAccelerator, mem = 15))
```