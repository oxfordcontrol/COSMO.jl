# # Custom Cone Constraint
#
# This example demonstrate how the user can extend `COSMO` with his own custom convex cones and use them in constraints. To make things easy, we will implement the simple cone of **Nonpositives**, i.e. $\mathbb{R}_{-}^n = \{x \in   \mathbb{R}^n \mid x_i \leq 0\}$.
# ## Cone definition and projection function
# We start by creating a concrete subtype of our cone interface-type: `COSMO.AbstractConvexCone{T}`, where `T` is a parameter for the floating-point precision. The only field that is strictly required for our cone definition is `dim`, which is used by the solver to store information about the cone dimension. However, you can also use the object to store other information, or allocate the workspace necessary for the projection step.
using COSMO, LinearAlgebra, SparseArrays, Test

## define new cone type
struct Nonpositives{T} <: COSMO.AbstractConvexCone{T}
    dim::Int64
end
#-
# Next, we define a projection method for `Nonpositives{T}` that takes a vector `x` and projects it onto our custom cone. For the cone of nonpositive vectors the projection simply sets all positive elements of `x` to zero:

function COSMO.project!(x::AbstractVector{T}, C::Nonpositives{T}) where {T <: AbstractFloat}
    x .= min.(x, zero(T))
end

#-
# This is all that is required to get a basic constraint working. Let's test our new cone/constraint by solving the following LP:
# $$
# \begin{array}{ll} \text{maximize} & x_1 + x_2 + x_3\\
# \text{subject to} &  x_1 \leq 3 \\
#                   &  x_2 \leq 2 \\
#                   &  x_1 + x_3 =  5
# \end{array}
# $$
# We define our decision vector $x = (x_1, x_2, x_3)$ and setup the problem:
n = 3
P = spzeros(n, n)
q = -ones(n)

## A1 * x + b1 <= 0
A1 = diagm(0 => ones(2))
b1 = [-3.; -2.]
c1 = COSMO.Constraint(A1, b1, Nonpositives, n, 1:2); #here we use our new cone

## x_1 + x_3 == 5
A2 = [1. 0 1.]
b2 = [-5.]
c2 = COSMO.Constraint(A2, b2, COSMO.ZeroSet);

## assemple and solve
model = COSMO.Model();
assemble!(model, P, q, [c1; c2]);
res = COSMO.optimize!(model);

#-
res.x

#-
obj_val = -dot(q, res.x)
# The problem was solved using constraints involving our new cone and yields the expected result. Next, we want to add the ability to detect infeasible problems.
#
# ## Support infeasibility detection
# If no further information about the new cone is provided, the infeasibility detection is disabled. However, by defining the two additional methods `in_dual` and `in_pol_recc` we can support infeasibility detection. More information on how infeasible problems are detected in `COSMO` can be found here \[1\]. We will have to add the ability for the solver to check whether a vector is in the [dual cone](https://en.wikipedia.org/wiki/Dual_cone_and_polar_cone) $\mathcal{K}^*$ of our new cone $\mathcal{K}$ and whether a vector is in the recession cone of the [polar cone](https://en.wikipedia.org/wiki/Dual_cone_and_polar_cone) ${\mathcal{K}^\circ}^\infty$ of our cone. The first function is used to check for primal infeasibility, while the second function  is used to check for dual infeasibility. Luckily, for the Nonpositives-cone this is straightforward:
#
# $$
# {\mathbb{R}_{-}^n}^* = \mathbb{R}_{-}^n, \quad {{\mathbb{R}_{-}^n}^\circ}^\infty = \mathbb{R}_{+}^n.
# $$

function COSMO.in_dual(x::AbstractVector{T}, C::Nonpositives{T}, tol::T) where {T <: AbstractFloat}
	return !any( x-> (x > -tol), x)
end

function COSMO.in_pol_recc(x::AbstractVector{T}, C::Nonpositives{T}, tol::T) where {T <: AbstractFloat}
	return !any( x-> (x < tol), x)
end
# Notice that for numerical reasons, we give the check-functions a bit of slack using the tolerance parameter `tol`.
# To test the infeasibility detection, we will attempt to use our new cone to solve the clearly dual infeasible problem:
# $$
# \begin{array}{ll} \text{minimize} & x\\
# \text{subject to} &  x\leq 3
# \end{array}
# $$
n = 1
P = spzeros(n, n)
q = [1.]

##  x <= 3 <=> x - 3 in Nonpositives
Ai = [1.]
bi = [-3.]
ci = COSMO.Constraint(Ai, bi, Nonpositives);

model = COSMO.Model();
assemble!(model, P, q, [ci]);
res = COSMO.optimize!(model);
#-
@test res.status == :Dual_infeasible

# which is what we would expect.
# ## Using new cone with JuMP
# The nice composibility of Julia-code allows us to use our new cone definition even when the problem is modelled with `JuMP`. For this to work, we define a concrete subtype of `MOI.AbstractSet`, that JuMP can use to build a constraint.
using MathOptInterface, JuMP
import Base.copy
const MOI = MathOptInterface

struct NonPos <: MOI.AbstractVectorSet
	dimension::Int
end
Base.copy(set::NonPos) = set #this is needed for MOI
# Next, we tell `COSMO` that whenever `JuMP` wants to solve a problem with a `NonPos`-cone it should translate it into a `Nonpositives` constraint and we also tell `MOI` that we now support constraints of this new mysterious cone `NonPos`.
function COSMO.processSet!(b::Vector{T}, rows::UnitRange{Int}, cs, s::NonPos) where {T <: AbstractFloat}
    push!(cs, Nonpositives{T}(length(rows)))
    nothing
end
#tell MOI, that COSMO supports constraints with this new cone
MOI.supports_constraint(optimizer::COSMO.Optimizer, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}}, ::Type{NonPos}) = true


# We should now be able to model the LP from above in `JuMP` and solve it internally using our `Nonpositives` cone and projection function.
model = JuMP.Model(COSMO.Optimizer);
@variable(model, x[1:3]);
@objective(model, Max, x[1] + x[2] + x[3]);
@constraint(model, A1 * x[1:2] .+ b1 in NonPos(2));
@constraint(model, x[1] + x[3] == 5);
JuMP.optimize!(model)

#-
x_opt = JuMP.value.(x)
# You can see in the solver output that indeed a set `Nonpositives` of `dim: 2` was used.
#
# The discussed cone of Nonpositives is of course trivial, but the ability to define new cones and constraints for `COSMO` can be very powerful to model complex problems.

# ## References
# [1] Garstka et al. - COSMO: A conic operator splitting method for convex conic problems (2020)
