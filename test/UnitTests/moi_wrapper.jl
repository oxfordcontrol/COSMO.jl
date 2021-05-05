using MathOptInterface, COSMO, Test, LinearAlgebra
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges
const MOIU = MOI.Utilities
MOIU.@model(COSMOModelData,
        (),
        (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
        (MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone,
         MOI.PositiveSemidefiniteConeTriangle, MOI.ExponentialCone, MOI.DualExponentialCone),
        (MOI.PowerCone, MOI.DualPowerCone),
        (),
        (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
        (MOI.VectorOfVariables,),
        (MOI.VectorAffineFunction,),);


struct UnsupportedModelAttribute  <: MOI.AbstractModelAttribute end

get_inner_optimizer(bridged) = bridged.model.optimizer

@testset "MOI Wrapper" begin

    @testset "Basic properties" begin
        optimizer =  COSMO.Optimizer(check_termination = 1, verbose = false);
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
    # mat(x) is posdeftriangle
    # with:
    A1_t = [1.0; 0; 3; 2; 14; 5];
    A2_t = [0.0; 4; 6; 16; 0; 4];
    C_t = [1.; 4; 9; 6; 0; 7];


    A1 = [1.0 0 1; 0 3 7; 1 7 5];
    A2 = [0.0 2 8; 2 6 0; 8 0 4];
    C = [1.0 2 3; 2 9 0; 3 0 7];
    b1 = 11.0;
    b2 = 19.0;

    
    cache = MOIU.UniversalFallback(MOIU.Model{Float64}());
    optimizer =  COSMO.Optimizer(check_termination = 1, verbose = false);
    cached = MOIU.CachingOptimizer(cache, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)

    x = MOI.add_variables(bridged, 6);
    # define objective function:
    objectiveFunction = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vec(C_t),x[1:6]),0.0);
    MOI.set(bridged, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),objectiveFunction);
    MOI.set(bridged, MOI.ObjectiveSense(), MOI.MIN_SENSE);

    # eq constraints
    con1 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A1_t,x[1:6]),0.0), MOI.EqualTo(b1));
    con2 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A2_t,x[1:6]),0.0), MOI.EqualTo(b2));

    # SDP constraint
    con3 = MOI.add_constraint(bridged, MOI.VectorOfVariables([x[1], x[2], x[3], x[4], x[5], x[6]]), MOI.PositiveSemidefiniteConeTriangle(3));

    # copy model into optimizer and solve
    @test MOI.get(bridged, MOI.ListOfVariableIndices()) == MOI.VariableIndex.(1:6)
    MOI.optimize!(bridged);
    @test sprint(show, optimizer) != nothing

    t_cold = MOI.get(bridged, MOI.SolveTime())
    iter_cold = MOI.get(bridged, COSMO.ADMMIterations())
    x_sol = MOI.get(bridged, MOI.VariablePrimal(), x)
    y_c1 = MOI.get(bridged, MOI.ConstraintDual(), con1)
    y_c2 = MOI.get(bridged, MOI.ConstraintDual(), con2)
    y_c3 = MOI.get(bridged, MOI.ConstraintDual(), con3)
    # Constraint Primal: For a constraint F-S, constraint primal is the value of the function evaluated at the primal solution
    s_c1 = MOI.get(bridged, MOI.ConstraintPrimal(), con1)
    s_c2 = MOI.get(bridged, MOI.ConstraintPrimal(), con2)
    s_c3 = MOI.get(bridged, MOI.ConstraintPrimal(), con3)

    @testset "SolverAttributes" begin
        @test MOI.get(bridged, MOI.SolverName()) == "COSMO";
        @test MOI.get(bridged, MOI.SolveTime()) >= 0.;
        # @test typeof(MOI.get(bridged, MOI.RawSolver())) <: COSMO.Workspace
        @test MOI.get(bridged, MOI.ResultCount()) == 1
        @test MOI.get(bridged, MOI.NumberOfVariables()) == 6
        @test isapprox(s_c1, b1, atol = 1e-3)
        @test isapprox(s_c2, b2, atol = 1e-3)
    end



    @testset "Warm starting" begin

        # solve the same problem again
        A1_t = [1.0; 0; 3; 2; 14; 5];
        A2_t = [0.0; 4; 6; 16; 0; 4];
        C_t = [1.; 4; 9; 6; 0; 7];


        cache = MOIU.UniversalFallback(MOIU.Model{Float64}());
        optimizer =  COSMO.Optimizer(check_termination = 1, verbose = false);
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, Float64)
        
        x = MOI.add_variables(bridged, 6);
        objectiveFunction = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(C_t, x[1:6]), 0.0);
        MOI.set(bridged, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), objectiveFunction);
        MOI.set(bridged, MOI.ObjectiveSense(), MOI.MIN_SENSE);
        con1 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A1_t, x[1:6]), 0.0), MOI.EqualTo(b1));
        con2 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A2_t, x[1:6]), 0.0), MOI.EqualTo(b2));
        con3 = MOI.add_constraint(bridged, MOI.VectorOfVariables([x[1], x[2], x[3], x[4], x[5], x[6]]), MOI.PositiveSemidefiniteConeTriangle(3));
        MOI.optimize!(bridged);

        x_sol_tri = MOI.get(bridged, MOI.VariablePrimal(), x)
        X_sol_tri = zeros(3, 3)
        optimizer = get_inner_optimizer(bridged)
        internal_scaled_s = copy(optimizer.inner.vars.s.data)
        internal_scaled_μ = copy(optimizer.inner.vars.μ)
        y_c3 = MOI.get(bridged, MOI.ConstraintDual(), con3)
        s_c3 = MOI.get(bridged, MOI.ConstraintPrimal(), con3)


        # warm start the PSD triangle related constraint primal and dual variables
        MOI.set(bridged, MOI.ConstraintDualStart(), con3, y_c3)
        MOI.set(bridged, MOI.ConstraintPrimalStart(), con3, s_c3)
        MOI.optimize!(bridged);
        
        optimizer = get_inner_optimizer(bridged) 
        # @test isapprox(optimizer.inner.vars.μ[3:end], internal_scaled_μ[3:end], atol = 1e-3)
        # @test isapprox(optimizer.inner.vars.s.data[3:end], internal_scaled_s[3:end], atol = 1e-4)

        # Warm start triggered by copy_to function (warm starting values from the model)
        # Different QP problem:
        # min 1/2 x'Px + q'x
        # with  P = [4. 1;1 2]; q = [1; 1.]
        # s.t. Ax <= u   <=> -Ax + u in Nonnegatives
        #      Ax >= l   <=> Ax - l in Nonnegatives
        # with A = [1. 1;1 0; 0 1];
        l = [1.; 0; 0];
        u = [1.; 0.7; 0.7]



        cache = MOIU.UniversalFallback(MOIU.Model{Float64}());
        optimizer =  COSMO.Optimizer();
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, Float64)
        MOI.set(optimizer, MOI.Silent(), true)


        x = MOI.add_variables(bridged, 2);
        objectiveFunction = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]); MOI.ScalarAffineTerm(1.0, x[2])], [MOI.ScalarQuadraticTerm(4.0, x[1], x[1]); MOI.ScalarQuadraticTerm(1.0, x[1], x[2]); MOI.ScalarQuadraticTerm(2.0, x[2], x[2])] , 0);
        MOI.set(bridged, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), objectiveFunction);
        MOI.set(bridged, MOI.ObjectiveSense(), MOI.MIN_SENSE);
        A = [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[2])),MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, x[1])),MOI.VectorAffineTerm(3, MOI.ScalarAffineTerm(1.0, x[2]))];
        con1 = MOI.add_constraint(bridged, MOI.VectorAffineFunction(A, -u), MOI.Nonpositives(3));
        con2 = MOI.add_constraint(bridged, MOI.VectorAffineFunction(A, -l),MOI.Nonnegatives(3));

        MOI.optimize!(bridged);

        # solve once to get optimal solution
        x_sol = MOI.get(bridged, MOI.VariablePrimal(), x);
        y_c1 = MOI.get(bridged, MOI.ConstraintDual(), con1);
        y_c2 = MOI.get(bridged, MOI.ConstraintDual(), con2);
        optimizer = get_inner_optimizer(bridged)
        y_c1_rows = optimizer.rowranges[1];
        y_c2_rows = optimizer.rowranges[2];

        # provide warm start values to the model
        MOI.set.(bridged, MOI.VariablePrimalStart(), x, x_sol)
        MOI.set.(bridged, MOI.ConstraintDualStart(), [con1, con2], [y_c1, y_c2])
        MOI.empty!(optimizer);


        MOIU.attach_optimizer(bridged.model)
;

        # check that internal variables are set correctly
        @test optimizer.inner.vars.x == x_sol
        @test optimizer.inner.vars.μ[y_c1_rows] == y_c1
        @test optimizer.inner.vars.μ[y_c2_rows] == -y_c2

    end

    @testset "Small edge cases" begin
        cache = MOIU.UniversalFallback(MOIU.Model{Float64}());
        optimizer =  COSMO.Optimizer(max_iter = 2, verbose = false);
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, Float64)
        
        x = MOI.add_variables(bridged, 1);
        objf = MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(-1.0, x[1])] , 0);
        MOI.set(bridged, MOI.ObjectiveFunction{typeof(objf)}(), objf);
        MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], x[1:1]), 0.0), MOI.GreaterThan{Float64}(10.));

        MOI.optimize!(bridged);
        @test MOI.get(bridged, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
        @test MOI.get(bridged, MOI.PrimalStatus()) == MOI.NEARLY_FEASIBLE_POINT
        @test MOI.get(bridged, MOI.DualStatus()) == MOI.NEARLY_FEASIBLE_POINT

    end
    @testset "Set merging" begin
        # based on HS21 from Maros Meszaros QP Set but with a range of redundant constraints to test set merging
        # the original problem is:
        # min 1/2 x' P x + r
        # s.t. l <= Ax <= u with the following data
        P =[0.02 0; 0 2.0]
        r = -100.;
        A = [-10.0 1; -1 0; 0 -1]
        l = [ 10.0; 2; -50];
        u = [1e20; 50.0; 50.0]

        # the known solution is
        obj_true = -99.96
        x_true = [-2.; 0]


    
        cache = MOIU.UniversalFallback(MOIU.Model{Float64}());
        optimizer =  COSMO.Optimizer(verbose = false);
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, Float64)

        x = MOI.add_variables(bridged, 2);
        objectiveFunction = MOI.ScalarQuadraticFunction{Float64}([], [MOI.ScalarQuadraticTerm(0.02, x[1], x[1]); MOI.ScalarQuadraticTerm(2.0, x[2], x[2])], r);
        MOI.set(bridged, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), objectiveFunction);
        MOI.set(bridged, MOI.ObjectiveSense(), MOI.MIN_SENSE);

        # Add the constraints in an unsorted order to make merging more "challenging"
        # A[1, :] * x >= 10
        nn1 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[1, :], x[1:2]), 0.0), MOI.GreaterThan(l[1]));
        # Redundant SOC constraint: -x[1] >= ||x[2]||_2
        soc = MOI.add_constraint(bridged, MOI.VectorAffineFunction(MOI.VectorAffineTerm.( [1; 2], MOI.ScalarAffineTerm.([-1. ;1. ], x[1:2])), zeros(2)), MOI.SecondOrderCone(2));
        # Interval: 2 <= A[2, :] * x <= 50
        box1 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[2, :], x[1:2]), 0.0), MOI.Interval(l[2], u[2]));
        # Equal To: x[1] == - 2
        zero1 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], x[1:1]), 0.0), MOI.EqualTo(x_true[1]));
        # Nonpositives: x[1; 2] - [10; 10] <= [0; 0]
        nn2 = MOI.add_constraint(bridged, MOI.VectorAffineFunction(MOI.VectorAffineTerm.( [1; 2], MOI.ScalarAffineTerm.([1. ;1. ], x[1:2])), -10. .* ones(2)), MOI.Nonpositives(2));
        # Interval: -50 <= A[3, :] * x <= 50
        box2 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[3, :], x[1:2]), 0.0), MOI.Interval(l[3], u[3]));
        # ZeroSet: x[1; 2] - x_true == [0; 0]
        zeroset = MOI.add_constraint(bridged, MOI.VectorAffineFunction(MOI.VectorAffineTerm.( [1; 2], MOI.ScalarAffineTerm.([1. ;1. ], x[1:2])), -x_true), MOI.Zeros(2));

        # copy model into optimizer
        MOI.optimize!(bridged);

        # check that the sets have been merged
        optimizer = get_inner_optimizer(bridged)
        C = optimizer.inner.p.C
        @test length(C.sets) == 4
        @test typeof(C.sets[1]) <: COSMO.ZeroSet
        @test typeof(C.sets[2]) <: COSMO.Nonnegatives
        @test typeof(C.sets[3]) <: COSMO.Box
        @test typeof(C.sets[4]) <: COSMO.SecondOrderCone
        x_sol = MOI.get(optimizer, MOI.VariablePrimal(), x)
        @test isapprox( MOI.get(optimizer, MOI.ObjectiveValue()), obj_true, atol = 1e-3)
        @test norm(x_sol - x_true, Inf) < 1e-3

     end

     @testset "Unsupported Attributes" begin

        cache = MOIU.UniversalFallback(MOIU.Model{Float64}());
        optimizer =  COSMO.Optimizer{Float64}();
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, Float64)
     
     

         x = MOI.add_variable(bridged)
         y = MOI.add_variable(bridged)

         cf = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm{Float64}.([1., 1.], [x, y]), 0.)
         c = MOI.add_constraint(bridged, cf, MOI.LessThan(0.))

         MOI.set(bridged, UnsupportedModelAttribute(), 0.)
         @test_throws MOI.UnsupportedAttribute{UnsupportedModelAttribute} MOI.optimize!(bridged);


     end

     @testset "BigFloat Test" begin
         T = BigFloat

         # min 1/2 x'Px + q'x
         # with  P = [4. 1;1 2]; q = [1; 1.]
         # s.t. Ax <= u   <=> -Ax + u in Nonnegatives
         #      Ax >= l   <=> Ax - l in Nonnegatives
         # with A = [1. 1;1 0; 0 1];
         l = T[1.; 0; 0];
         u = T[1.; 0.7; 0.7]


        cache = MOIU.UniversalFallback(MOIU.Model{T}());
        optimizer =  COSMO.Optimizer{T}(verbose = false);
        cached = MOIU.CachingOptimizer(cache, optimizer)
        bridged = MOIB.full_bridge_optimizer(cached, T)

        MOI.set(bridged, MOI.Silent(), true)

         x = MOI.add_variables(bridged, 2);
         objectiveFunction = MOI.ScalarQuadraticFunction{T}([MOI.ScalarAffineTerm(T(1.0), x[1]); MOI.ScalarAffineTerm(T(1.0), x[2])], [MOI.ScalarQuadraticTerm(T(4.0), x[1], x[1]); MOI.ScalarQuadraticTerm(T(1.0), x[1], x[2]); MOI.ScalarQuadraticTerm(T(2.0), x[2], x[2])] , zero(T));
         MOI.set(bridged, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(), objectiveFunction);
         MOI.set(bridged, MOI.ObjectiveSense(), MOI.MIN_SENSE);
         A = [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(T(1.0), x[1])),MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(T(1.0), x[2])),MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(T(1.0), x[1])),MOI.VectorAffineTerm(3, MOI.ScalarAffineTerm(T(1.0), x[2]))];
         con1 = MOI.add_constraint(bridged, MOI.VectorAffineFunction(A, -u), MOI.Nonpositives(3));
         con2 = MOI.add_constraint(bridged, MOI.VectorAffineFunction(A, -l),MOI.Nonnegatives(3));

         MOI.optimize!(bridged);
         optimizer = get_inner_optimizer(bridged)
         @test optimizer.results.status == :Solved
         @test typeof(optimizer.results) == COSMO.Result{BigFloat}
         @test isapprox(optimizer.results.obj_val, 1.88, atol = 1e-3)
     end


    # # ------------------
    # # MOI - Test sets
    # --------------------

    optimizer = COSMO.Optimizer(eps_abs = 1e-5, eps_rel = 1e-5, check_termination = 25, check_infeasibility = 50)
    MOI.set(optimizer, MOI.Silent(), true)

    cache = MOIU.UniversalFallback(MOIU.Model{Float64}())
    cached = MOIU.CachingOptimizer(cache, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)
    config = MOIT.TestConfig(atol = 1e-3, rtol = 1e-3, duals = false)


    @testset "Unit" begin
        MOIT.unittest(bridged, config,
                        [# Quadratic constraints are not supported
                        "solve_qcp_edge_cases",
                       # ArgumentError: The number of constraints in SCSModel must be greater than 0
                        "solve_unbounded_model",
                       # Integer and ZeroOne sets are not supported
                       "solve_integer_edge_cases", "solve_objbound_edge_cases",
                        "solve_zero_one_with_bounds_1",
                        "solve_zero_one_with_bounds_2",
                        "solve_zero_one_with_bounds_3",
                       # Can not set the # of threads
                        "number_threads",
                        "solve_farkas_interval_upper"
                         ])
    end

    @testset "Continuous Linear" begin
        MOIT.contlineartest(bridged, config, [
            "linear8a", # It expects `ResultCount` to be 0 as we disable `duals`.
            ])
    end

    @testset "Continuous Quadratic" begin

        MOIT.contquadratictest(bridged, config, [ "socp",
                                                 "ncqcp"])
    end


    exclude_conic_test_sets = ["rootdet", "logdet"]
    @testset "Continuous Conic" begin
        MOIT.contconictest(bridged, config, exclude_conic_test_sets)
    end

 end
