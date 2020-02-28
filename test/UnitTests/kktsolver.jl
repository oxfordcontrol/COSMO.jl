using COSMO, Test, LinearAlgebra, SparseArrays, Random, QDLDL, Pkg
include("COSMOTestUtils.jl")

rng = Random.MersenneTwister(1)

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

@testset "KKT Solver" begin
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
        if isdefined(Pardiso, :MKL_PARDISO_LOADED)
            push!(solver_types, COSMO.MKLPardisoKKTSolver)
            push!(solver_tols, 1e-10)
                add_kwargs(params)
            end
        if isdefined(Pardiso, :PARDISO_LOADED)
            push!(solver_types, COSMO.PardisoDirectKKTSolver)
            push!(solver_tols, 1e-10)
            add_kwargs(params)
        end
        if isdefined(Pardiso, :PARDISO_LOADED)
            push!(solver_types, COSMO.PardisoIndirectKKTSolver)
            push!(solver_tols, 5e-5)
            add_kwargs(params)
        end
    end

    if in("IterativeSolvers",keys(Pkg.installed())) && in("LinearMaps",keys(Pkg.installed()))
        using IterativeSolvers
        push!(solver_types, COSMO.IndirectReducedKKTSolver)
        push!(solver_tols, 1e-3)
        add_kwargs(params, solver_type=:CG)

        push!(solver_types, COSMO.IndirectReducedKKTSolver)
        push!(solver_tols, 1e-3)
        add_kwargs(params, solver_type=:MINRES)

        push!(solver_types, COSMO.IndirectKKTSolver)
        push!(solver_tols, 1e-3)
        add_kwargs(params, solver_type=:MINRES)
    end


    @testset "$(solver_types[i]): KKT solver tests" for i = 1:length(solver_types)

        m = 10
        n = 20

        for rho1 in [rand(T), rand(T, m)],
            rho2 in [rand(T), rand(T, m)],
            sigma in [rand(T)]

            P  = sparse(generate_pos_def_matrix(rng, n))
            A  = sprandn(m, n, 0.2)
            b = randn(m + n)

            F = solver_types[i](P, A, sigma, rho1; params[i]...)
            if in("IterativeSolvers",keys(Pkg.installed())) && in("LinearMaps",keys(Pkg.installed()))
                if isa(F, COSMO.IndirectReducedKKTSolver) || isa(F, COSMO.IndirectKKTSolver)
                    # Technically, we should have been able to set the tolerance even lower
                    # But, currently, issues in IterativeSolvers.jl do not allow this
                    # See: https://github.com/JuliaMath/IterativeSolvers.jl/pull/244
                    F.iteration_counter = 10^4 # This forces CG's/MINRES tolerance to be 1/10^4
                end
            end
            J = make_test_kkt(P, A, sigma, rho1)
            x = copy(b)

            #test a single solve
            COSMO.solve!(F, x, b)
            @test norm(x - J \ b) <= solver_tols[i]

            # Check that warm starting works
            # Invoking again an indirect solver should result in the solution with only
            # one matrix vector multiplication
            if in("IterativeSolvers",keys(Pkg.installed())) && in("LinearMaps",keys(Pkg.installed()))
                if isa(F, COSMO.IndirectReducedKKTSolver) || isa(F, COSMO.IndirectKKTSolver)
                    # The calculation of the residual, and thus the termination criterion, of
                    # MINRES is approximate. Thus warm started solutions won't necessarily finish in one step
                    # For this reason we don't check warm starting with MINRES for now :(
                    if F.solver_type != :MINRES
                        COSMO.solve!(F, x, b)
                        @test F.multiplications[end] <= 1
                        @test norm(x - J \ b) <= solver_tols[i]
                    end
                end
            end

            #test a rho update and solve
            J = make_test_kkt(P, A, sigma, rho2)
            x = copy(b)
            COSMO.update_rho!(F, rho2)
            COSMO.solve!(F, x, b)
            @test norm(x - J \ b) <= solver_tols[i]

         end

    end

    # Try every solver on the same problem and compare results
    solver_types = Array{Union{Type{<: COSMO.AbstractKKTSolver}, COSMO.OptionsFactory{<: COSMO.AbstractKKTSolver}}}(undef, 0)
    push!(solver_types, COSMO.QdldlKKTSolver, COSMO.CholmodKKTSolver)
    if in("Pardiso",keys(Pkg.installed()))
        isdefined(Pardiso, :MKL_PARDISO_LOADED) && push!(solver_types, COSMO.PardisoDirectKKTSolver)
        if isdefined(Pardiso, :PARDISO_LOADED)
            push!(solver_types, COSMO.PardisoDirectKKTSolver)
            push!(solver_types, COSMO.PardisoIndirectKKTSolver)
        end
    end
    if in("IterativeSolvers",keys(Pkg.installed())) && in("LinearMaps",keys(Pkg.installed()))
        push!(solver_types, COSMO.CGIndirectKKTSolver)
        push!(solver_types, COSMO.MINRESIndirectKKTSolver)
        push!(solver_types, with_options(COSMO.CGIndirectKKTSolver))
        push!(solver_types, with_options(COSMO.MINRESIndirectKKTSolver))
    end

     P = [4. 1;1 2]
     q = [1; 1.]
     A = [1. 1;1 0; 0 1]
     l = [1.;0;0]
     u = [1.;0.7;0.7]

     constraint1 = COSMO.Constraint(-A, u, COSMO.Nonnegatives)
     constraint2 = COSMO.Constraint(A, -l, COSMO.Nonnegatives)
     constraints = [constraint1; constraint2]
     @testset "$(solver_types[i]): KKT solver tests - simple problem" for i = 1:length(solver_types)
         model = COSMO.Model()
         settings = COSMO.Settings(kkt_solver = solver_types[i])
         assemble!(model, P, q, constraints, settings = settings)
         res = COSMO.optimize!(model);
         @test isapprox(norm(res.x - [0.3; 0.7]), 0., atol=1e-3)
         @test isapprox(res.obj_val, 1.88, atol=1e-3)

     end

end

nothing
