# Example from JuliaOpt/SumOfSquares.jl
# COSMO supports MOI/JuMP and therefore also works with JuMP extensions
# like SumOfSquares.jl



using SumOfSquares, COSMO
using DynamicPolynomials
@polyvar x y
motzkin = x^4*y^2 + x^2*y^4 + 1 - 3x^2*y^2

model = SOSModel(with_optimizer(COSMO.Optimizer, rho = 1e-5);)
@constraint(model, (x^2 + y^2) * motzkin >= 0)
optimize!(model)

primal_status(model)

# The full example can be found at:
# https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/motzkin.ipynb


