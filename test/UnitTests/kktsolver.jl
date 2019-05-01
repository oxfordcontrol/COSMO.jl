using COSMO, Test, LinearAlgebra, SparseArrays, Random, QDLDL, Pkg
include("COSMOTestUtils.jl")

rng = Random.MersenneTwister(2401)

T = Float64

function make_test_kkt(P, A, sigma, rho)

    R = length(rho)   == 1 ? ((1.) ./ rho[1]) * I : Diagonal((1.) ./ rho)
    S = length(sigma) == 1 ? (sigma[1]) * I : Diagonal(sigma)

    #compute the full KKT matrix
    K    = [P+S A'; A -R]
    return K
end

function add_kwargs(array; kwargs...)
    push!(array, kwargs)
end

solver_types = []
solver_tols = []
params = []
push!(solver_types, COSMO.QdldlKKTSolver)
push!(solver_tols, 1e-10)
add_kwargs(params)
push!(solver_types, COSMO.CholmodKKTSolver)
push!(solver_tols, 1e-10)
add_kwargs(params)
# solver_types = [COSMO.QdldlKKTSolver
                # COSMO.CholmodKKTSolver]

# optional dependencies
if in("Pardiso",keys(Pkg.installed()))
    using Pardiso
    if Pardiso.MKL_PARDISO_LOADED[] 
        push!(solver_types, COSMO.MKLPardisoKKTSolver)
        push!(solver_tols, 1e-10)
            add_kwargs(params)
        end
    if Pardiso.PARDISO_LOADED[] 
        push!(solver_types, COSMO.PardisoDirectKKTSolver)
        push!(solver_tols, 1e-10)
        add_kwargs(params)
    end
    if Pardiso.PARDISO_LOADED[]
        push!(solver_types, COSMO.PardisoIndirectKKTSolver)
        push!(solver_tols, 5e-5)
        add_kwargs(params)
    end
end

if in("IterativeSolvers",keys(Pkg.installed()))
    using IterativeSolvers
    push!(solver_types, COSMO.IndirectReducedKKTSolver)#{SparseMatrixCSC{T, Int}, SparseMatrixCSC{T, Int}, T})
    push!(solver_tols, 1e-6)
    add_kwargs(params, solver=:CG)

    push!(solver_types, COSMO.IndirectReducedKKTSolver)#{SparseMatrixCSC{T, Int}, SparseMatrixCSC{T, Int}, T})
    push!(solver_tols, 1e-6)
    add_kwargs(params, solver=:MINRES)
end


 @testset "$(solver_types[i]) $(params[i]) : KKT solver tests" for i = 1:length(solver_types)

    m = 10
    n = 20

    for rho1 in [rand(T), rand(T, m)],
        rho2 in [rand(T), rand(T, m)],
        sigma in [rand(T)]

        P  = sparse(generate_pos_def_matrix(n, rng))
        A  = sprandn(m, n, 0.2)
        b = randn(m + n)

        F = solver_types[i](P, A, sigma, rho1; params[i]...)
        if isa(F, COSMO.IndirectReducedKKTSolver)
            F.iteration_counter = 10^8 # This forces the solution's tolerance to be very small
        end
        J = make_test_kkt(P, A, sigma, rho1)
        x = copy(b)

        #test a single solve
        COSMO.solve!(F, x, b)
        @test norm(x - J \ b) <= solver_tols[i]

        #test a rho update and solve
        J = make_test_kkt(P, A, sigma, rho2)
        x = copy(b)
        COSMO.update_rho!(F, rho2)
        COSMO.solve!(F, x, b)
        @test norm(x - J \ b) <= solver_tols[i]

     end

 end

nothing
