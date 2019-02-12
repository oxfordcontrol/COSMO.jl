using MathOptInterface, COSMO, Test, LinearAlgebra
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model(COSMOModelData,
        (),
        (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
        (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
         MOI.PositiveSemidefiniteConeSquare, MOI.PositiveSemidefiniteConeTriangle),
        (),
        (MOI.SingleVariable,),
        (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
        (MOI.VectorOfVariables,),
        (MOI.VectorAffineFunction,),);


@testset "MOI Wrapper" begin

    @testset "Basic properties" begin
        optimizer =  COSMO.Optimizer(check_termination = 1);
        @test sprint(show, optimizer) == string("Empty COSMO - Optimizer")
        @test MOI.supports(optimizer, MOI.ObjectiveFunction{MOI.SingleVariable}())
        @test MOI.supports(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
        @test MOI.supports(optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}())
        @test MOI.supports(optimizer, MOI.ObjectiveSense())
        @test MOI.supports(optimizer, MOI.VariablePrimalStart(), MOI.VariableIndex)
        @test MOI.supports(optimizer, MOI.ConstraintPrimalStart(), MOI.ConstraintIndex)
        @test MOI.supports(optimizer, MOI.ConstraintDualStart(), MOI.ConstraintIndex)
    end
    # Solve the following problem:
    # min c'*x
    # vec(A1)' * x == b1
    # vec(A2)' * x == b2
    # mat(x) is posdef
    # with:
    A1 = [1.0 0 1; 0 3 7; 1 7 5];
    A2 = [0.0 2 8; 2 6 0; 8 0 4];
    C = [1.0 2 3; 2 9 0; 3 0 7];
    b1 = 11.0;
    b2 = 19.0;

    model = COSMOModelData{Float64}();
    optimizer =  COSMO.Optimizer(check_termination = 1);
    x = MOI.add_variables(model, 9);
    # define objective function:
    objectiveFunction = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vec(C),x[1:9]),0.0);
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),objectiveFunction);
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE);

    # eq constraints
    con1 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vec(A1),x[1:9]),0.0), MOI.EqualTo(b1));
    con2 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vec(A2),x[1:9]),0.0), MOI.EqualTo(b2));

    # SDP constraint
    con3 = MOI.add_constraint(model, MOI.VectorOfVariables(x[1:9]), MOI.PositiveSemidefiniteConeSquare(3));

    # copy model into optimizer
    MOI.empty!(optimizer);
    @test sprint(show, optimizer) != nothing

    idxmap = MOI.copy_to(optimizer, model);
    @test MOI.get(optimizer, MOI.ListOfVariableIndices()) == MOI.VariableIndex.(1:9)
    MOI.optimize!(optimizer);
    @test sprint(show, optimizer) != nothing

    t_cold = MOI.get(optimizer, MOI.SolveTime())
    iter_cold = optimizer.results.iter
    x_sol = MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), x))
    X_sol = reshape(x_sol, 3, 3)
    y_c1 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con1])
    y_c2 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con2])
    y_c3 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con3])
    # Constraint Primal: For a constraint F-S, constraint primal is the value of the function evaluated at the primal solution
    s_c1 = MOI.get(optimizer, MOI.ConstraintPrimal(), idxmap[con1])
    s_c2 = MOI.get(optimizer, MOI.ConstraintPrimal(), idxmap[con2])
    s_c3 = MOI.get(optimizer, MOI.ConstraintPrimal(), idxmap[con3])

    @testset "SolverAttributes" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "COSMO";
        @test MOI.get(optimizer, MOI.SolveTime()) > 0.;
        @test typeof(MOI.get(optimizer, MOI.RawSolver())) <: COSMO.Workspace
        @test MOI.get(optimizer, MOI.ResultCount()) == 1
        @test MOI.get(optimizer, MOI.NumberOfVariables()) == 9
        @test isapprox(s_c1, b1, atol = 1e-3)
        @test isapprox(s_c2, b2, atol = 1e-3)
        # check if S_C3 = X is pos semidefinite
        @test minimum(eigen(reshape(s_c3, 3, 3)).values) >  -1e-6
        # check if Y_C3 is pos semidefinite (since y is in dual cone K*)
        @test minimum(eigen(reshape(y_c3, 3, 3)).values) >  -1e-6
    end

    # Solve once again to get cold iter and time
    MOI.empty!(optimizer);
    idxmap = MOI.copy_to(optimizer, model);
    MOI.optimize!(optimizer);
    iter_cold = optimizer.results.iter

    @testset "Warm starting" begin

        # Warm start x,y, s at the same time
        optimizer =  COSMO.Optimizer(check_termination = 1);
        MOI.empty!(optimizer);
        copyresult = MOI.copy_to(optimizer, model);
        MOI.set.(optimizer, MOI.VariablePrimalStart(), x[1:9], x_sol[1:9])
        MOI.set.(optimizer, MOI.ConstraintPrimalStart(), [con1, con2, con3], [s_c1, s_c2, s_c3])
        MOI.set.(optimizer, MOI.ConstraintDualStart(), [con1, con2, con3], [y_c1, y_c2, y_c3])

        # check that variables are correctly set after the warm starting
        @test x_sol == optimizer.inner.vars.x
        @test y_c1 == optimizer.inner.vars.μ[1]
        @test y_c2 == optimizer.inner.vars.μ[2]
        @test y_c3 == -optimizer.inner.vars.μ[3:end]
        @test 0. == optimizer.inner.vars.s.data[1]
        @test 0. == optimizer.inner.vars.s.data[2]
        @test s_c3 == optimizer.inner.vars.s.data[3:end]

        MOI.optimize!(optimizer);
        iter_warm_all = optimizer.results.iter
        @test iter_warm_all < iter_cold

        # solve the same problem but with a psd triangle constraint (upper triangle)
        A1_t = [1.0; 0; 3; 2; 14; 5];
        A2_t = [0.0; 4; 6; 16; 0; 4];
        C_t = [1.; 4; 9; 6; 0; 7];
        model = COSMOModelData{Float64}();
        x = MOI.add_variables(model, 6);
        objectiveFunction = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(C_t, x[1:6]), 0.0);
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), objectiveFunction);
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE);
        con1 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A1_t, x[1:6]), 0.0), MOI.EqualTo(b1));
        con2 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A2_t, x[1:6]), 0.0), MOI.EqualTo(b2));
        con3 = MOI.add_constraint(model, MOI.VectorOfVariables([x[1], x[2], x[3], x[4], x[5], x[6]]), MOI.PositiveSemidefiniteConeTriangle(3));
        MOI.empty!(optimizer);

        idxmap = MOI.copy_to(optimizer, model);
        MOI.optimize!(optimizer);

        # double check solution with problem with PSDSquare constraint
        x_sol_tri = MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), x))
        X_sol_tri = zeros(3, 3)
        internal_scaled_s = copy(optimizer.inner.vars.s.data)
        internal_scaled_μ = copy(optimizer.inner.vars.μ)
        COSMO.populate_upper_triangle!(X_sol_tri, x_sol_tri, 1.)
        @test maximum(abs.(X_sol - Symmetric(X_sol_tri))) < 1e-2

        y_c3 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con3])
        s_c3 = MOI.get(optimizer, MOI.ConstraintPrimal(), idxmap[con3])

        MOI.empty!(optimizer);
        idxmap = MOI.copy_to(optimizer, model);

        # warm start the PSD triangle related constraint primal and dual variables
        MOI.set(optimizer, MOI.ConstraintDualStart(), con3, y_c3)
        MOI.set(optimizer, MOI.ConstraintPrimalStart(), con3, s_c3)
        @test isapprox(optimizer.inner.vars.μ[3:end], internal_scaled_μ[3:end], atol = 1e-6)
        @test isapprox(optimizer.inner.vars.s.data[3:end], internal_scaled_s[3:end], atol = 1e-6)

        # Warm start triggered by copy_to function (warm starting values from the model)
        # Different QP problem:
        # min 1/2 x'Px + q'x
        # with  P = [4. 1;1 2]; q = [1; 1.]
        # s.t. Ax <= u   <=> -Ax + u in Nonnegatives
        #      Ax >= l   <=> Ax - l in Nonnegatives
        # with A = [1. 1;1 0; 0 1];
        l = [1.; 0; 0];
        u = [1.; 0.7; 0.7]

        model = MOIU.UniversalFallback(COSMOModelData{Float64}());
        optimizer =  COSMO.Optimizer();

        x = MOI.add_variables(model, 2);
        objectiveFunction = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]); MOI.ScalarAffineTerm(1.0, x[2])], [MOI.ScalarQuadraticTerm(4.0, x[1], x[1]); MOI.ScalarQuadraticTerm(1.0, x[1], x[2]); MOI.ScalarQuadraticTerm(2.0, x[2], x[2])] , 0);
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), objectiveFunction);
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE);
        A = [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[2])),MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, x[1])),MOI.VectorAffineTerm(3, MOI.ScalarAffineTerm(1.0, x[2]))];
        con1 = MOI.add_constraint(model, MOI.VectorAffineFunction(A, -u), MOI.Nonpositives(3));
        con2 = MOI.add_constraint(model, MOI.VectorAffineFunction(A, -l),MOI.Nonnegatives(3));

        MOI.empty!(optimizer);
        idxmap = MOI.copy_to(optimizer, model);
        MOI.optimize!(optimizer);

        # solve once to get optimal solution
        x_sol = MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), x));
        y_c1 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con1]);
        y_c2 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con2]);
        y_c1_rows = COSMO.constraint_rows(optimizer.rowranges, idxmap[con1]);
        y_c2_rows = COSMO.constraint_rows(optimizer.rowranges, idxmap[con2]);

        # provide warm start values to the model
        MOI.set.(model, MOI.VariablePrimalStart(), x, x_sol)
        MOI.set.(model, MOI.ConstraintDualStart(), [con1, con2], [y_c1, y_c2])
        MOI.empty!(optimizer);
        idxmap = MOI.copy_to(optimizer, model);

        # check that internal variables are set correctly
        @test optimizer.inner.vars.x == x_sol
        @test optimizer.inner.vars.μ[y_c1_rows] == y_c1
        @test optimizer.inner.vars.μ[y_c2_rows] == -y_c2

    end

    @testset "Small edge cases" begin
        model = COSMOModelData{Float64}();
        optimizer =  COSMO.Optimizer(max_iter = 2);

        x = MOI.add_variables(model, 1);
        objectiveFunction = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(-1.0, x[1])] , 0);
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), objectiveFunction);
        MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.GreaterThan{Float64}(10.));
        idxmap = MOI.copy_to(optimizer, model);
        MOI.optimize!(optimizer);
        @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
        @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.NO_SOLUTION
        @test MOI.get(optimizer, MOI.DualStatus()) == MOI.NO_SOLUTION

    end
end
# --------------------
# MOI - Test sets
# --------------------
optimizer = MOIU.CachingOptimizer(MOIU.UniversalFallback(COSMOModelData{Float64}()),
                                  COSMO.Optimizer(eps_abs = 1e-6, eps_rel = 1e-6 ));


config = MOIT.TestConfig(atol=1e-2, rtol=1e-2)
@testset "Unit" begin
    MOIT.unittest(MOIB.SplitInterval{Float64}(optimizer), config,
                  [# Quadratic functions are not supported
                   "solve_qcp_edge_cases", "solve_qp_edge_cases",
                   # Integer and ZeroOne sets are not supported
                   "solve_integer_edge_cases", "solve_objbound_edge_cases"])
end

# FIXME: one test fails for linear12 problem. Take a closer look why
@testset "Continuous linear problems" begin
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config, ["linear12"])
end

@testset "Continuous quadratic problems" begin
    exclude_qt_test_sets = ["qcp", "socp"]
    MOIT.contquadratictest(optimizer, config, exclude_qt_test_sets)
end


exclude_conic_test_sets = ["exp", "rootdet", "logdet"]
@testset "Continuous conic problems" begin
    MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.LogDet{Float64}(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(optimizer)))),
                       config, exclude_conic_test_sets)
end
nothing