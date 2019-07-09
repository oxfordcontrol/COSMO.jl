# Closest Correlation matrix test matrix


# Original problem has the following format:
# min_X   1/2 ||X-C||^2
# s.t.    Xii = 1
#         X ⪴ 0

include("../../src/COSMO.jl")
using Main.COSMO, Test, LinearAlgebra, SparseArrays, Random


# this function creates a matrix A that slices out the diagonal entries Xii of a vectorized square matrix x=vec(X)
function create_diagonal_extractor(n)
    A = spzeros(n, Int(n*(n + 1)/2))
    idx = 0
    for i = 1:n
        idx += i
        A[i, idx] = 1
    end
    return A
end

function diagonal_equal_one(A, atol)
    for iii = 1:size(A, 1)
        (abs(A[iii, iii] - 1) > atol) && return false
    end
    return true
end

rng = Random.MersenneTwister(12345)
xMin = -1.
xMax = 1.
n = 200 
C = xMin .+ randn(rng, n, n)*(xMax - xMin)
C = Symmetric(C)

isposdef(C) && warn("The perturbed correlation matrix is still pos def.")

n2 = Int(n * (n + 1)/2)
m = n + n2
P = sparse(1.0I, n2, n2)
q = -COSMO.extract_upper_triangle(C)
r = 0.5*dot(q, q)

A1 = create_diagonal_extractor(n)
b1 = -ones(n)
A2 = sparse(1.0I, n2, n2)
b2 = zeros(n2)

function solve_problem(C, n, lanczos=true)
    cs1 = COSMO.Constraint(A1, b1, COSMO.ZeroSet)
    if lanczos
        cs2 = COSMO.Constraint(A2, b2, COSMO.PsdConeTriangleLanczos)
    else
        cs2 = COSMO.Constraint(A2, b2, COSMO.PsdConeTriangle)
    end
    constraints = [cs1; cs2]


    settings = COSMO.Settings(verbose=true, verbose_timing=true)
    model = COSMO.Model()
    COSMO.assemble!(model, P, q, constraints, settings=settings)

    res = COSMO.optimize!(model)
    return res, model
end

# @show Xsol

#=
using Profile, ProfileView
Profile.init(n = 10^7, delay = 0.0001)
Profile.clear()		
Profile.@profile res, model = solve_problem(C, n)
Profile.clear()		
Profile.@profile res, model = solve_problem(C, n)
Profile.clear()		
Profile.@profile res, model = solve_problem(C, n)
ProfileView.view()		
println("Press enter to continue...")		
readline(stdin)	
=#

using BenchmarkTools

res, model = solve_problem(C, n)
res, model = solve_problem(C, n)
solve_problem(C, n, false)
solve_problem(C, n, false)

Xsol = COSMO.populate_upper_triangle(res.x)
set = model.p.C.sets[2]
@show set.subspace_dim_history
@testset "Closest Correlation Matrix - SDP Problems" begin
@test res.status == :Solved
@test diagonal_equal_one(Xsol, 1e-5)
@test minimum(eigen(Xsol).values) > -1e-3

end
nothing
