# Optimal Eigenvalue problem described in Pataki 1998:
# On the Rank of Extreme Matrices in Semidefinite Programs and the Multiplicity of Optimal Eigenvalues

using COSMO, Test, LinearAlgebra, Random, SparseArrays
using JuMP, SCS, LinearAlgebra

n = 700;
m = 50;
k = 10;
rng = Random.MersenneTwister(12345)
A = []
for i = 1:m+1
    A_i = sprandn(rng, n, n, 0.05/2); A_i = (A_i + A_i')/2
    push!(A, A_i)
end
#@testset "Lanczos Projection - Simple SDP" begin
eye = Matrix{Float64}(I, n, n)
c = [eye[:]; zeros(n^2); k; zeros(m)]; c = reshape(c, length(c), 1)
P = spzeros(length(c), length(c))
speye = spdiagm(0 => ones(n^2))
speye_small = spdiagm(0 => ones(n))
global A_1 = [speye -speye speye_small[:]]
for i = 1:m
    global A_1 = [A_1 -A[i+1][:]]
end
b_1 = -Array(A[1][:]); b_1 = reshape(b_1, length(b_1), 1)
zero_constraint = COSMO.Constraint(A_1, b_1, COSMO.ZeroSet)
A_ = [speye spzeros(n^2, n^2 + m + 1)]
b_ = zeros(size(A_, 1))
conic_constraint_1 = COSMO.Constraint(A_, b_, COSMO.PsdCone)
A_ = [spzeros(n^2, n^2) speye spzeros(n^2, m + 1)]
b_ = zeros(size(A_, 1), 1)
conic_constraint_2 = COSMO.Constraint(A_, b_, COSMO.PsdCone)

# define example problem
settings = COSMO.Settings(verbose=true, verbose_timing=true, max_iter=2500)

model = COSMO.Model()
assemble!(model,P,c,[zero_constraint; conic_constraint_1; conic_constraint_2])

# res = COSMO.optimize!(model,settings)
# res = COSMO.optimize!(model,settings)

COSMO.optimize!(model,settings)
# @time COSMO.optimize!(model,settings)
using Profile
Profile.clear()
Profile.@profile COSMO.optimize!(model,settings)
using ProfileView
ProfileView.view()
@show res.times
println("Press enter to continue...")
readline(stdin)
#=
open("prof.txt", "w") do s
	Profile.print(IOContext(s, :displaysize => (24, 500)))
end
=#


model = JuMP.Model()
setsolver(model, SCSSolver(max_iters=10000, acceleration_lookback=1))
@variable(model, V[1:n,1:n], SDP)
@variable(model, W[1:n,1:n], SDP)
@variable(model, z)
@variable(model, x[1:m])
@objective(model, Min, k*z + sum(diag(V)))
@constraint(model, z.*speye_small + V - W .== A[1] + x[1].*A[2] + x[2].*A[3] + x[3].*A[4] + x[4].*A[5] + x[5].*A[6])
status = JuMP.solve(model)

println("The following two should be approximately the same")
@show getvalue(z)
@show res.x[2*n^2+1]

nothing