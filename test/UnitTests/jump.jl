# JuMP unit tests for MOI Wrapper
using COSMO, Test, LinearAlgebra, JuMP, MathOptInterface

# ---------------------------
# Linear Program
# ---------------------------

c = [300.0; 500.0]

m = Model(with_optimizer(COSMO.Optimizer))
@variable(m, x1 >= 0)
@variable(m, x2 >= 0)

@objective(m, Max, c[1] * x1 + c[2] * x2)
@constraint(m, x1 + 2 * x2 <= 170)
@constraint(m, x1 + x2 <= 150)
@constraint(m, 3 * x2 <= 180)
status = JuMP.optimize!(m)


termination_status = JuMP.termination_status(m)
result_status = JuMP.primal_status(m)
optimal_x1 = JuMP.value(x1)
optimal_x2 = JuMP.value(x2)
objvalue = JuMP.objective_value(m)

@testset "LP - JuMP" begin
    @test isapprox(JuMP.objective_value(m), 49000, atol = 1)
    @test isapprox(JuMP.value(x1), 130, atol = 0.1)
    @test isapprox(JuMP.value(x2), 20, atol = 0.1)
    @test JuMP.termination_status(m) == MathOptInterface.Success
    @test JuMP.primal_status(m) == MathOptInterface.FeasiblePoint
end


# ---------------------------
# Quadratic Program
# ---------------------------

m = Model(with_optimizer(COSMO.Optimizer))
@variable(m, x[1:2])

P = [4. 1; 1 2]
q = [1; 1.]
@objective(m, Min, 0.5 * x' * P * x  + q' * x)
A = [1. 1; 1 0; 0 1];
l = [1.; 0; 0];
u = [1.; 0.7; 0.7]
@constraint(m, A * x .<= u)
@constraint(m, A * x .>= l)
status = JuMP.optimize!(m)

@testset "QP - JuMP" begin
    @test isapprox(JuMP.objective_value(m), 1.88, atol = 1e-3)
    @test isapprox(JuMP.value.(x)[1], 0.3, atol = 1e-3)
    @test isapprox(JuMP.value.(x)[2], 0.7, atol = 1e-3)
    @test JuMP.termination_status(m) == MathOptInterface.Success
    @test JuMP.primal_status(m) == MathOptInterface.FeasiblePoint
end

# ---------------------------
# Second order cone Program
# ---------------------------

m = Model(with_optimizer(COSMO.Optimizer));
@variable(m, x[1:3]);
@objective(m, Min, x[1]);
@constraint(m, x[2] >= 1);
@constraint(m, x[3] >= 1);
@constraint(m,  soc, [x[1];x[2:3]] in SecondOrderCone());
status = JuMP.optimize!(m);

@testset "SOCP - JuMP" begin
    @test isapprox(JuMP.objective_value(m), sqrt(2), atol = 1e-3)
    @test isapprox(JuMP.value.(x)[1], sqrt(2), atol = 1e-3)
    @test isapprox(JuMP.value.(x)[2], 1., atol = 1e-3)
    @test isapprox(JuMP.value.(x)[3], 1., atol = 1e-3)
    @test JuMP.termination_status(m) == MathOptInterface.Success
    @test JuMP.primal_status(m) == MathOptInterface.FeasiblePoint
end

# ---------------------------
# Semidefinite Program - Square
# ---------------------------
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

m = Model(with_optimizer(COSMO.Optimizer));
@variable(m, X[1:3,1:3]);
@objective(m, Min, tr(C * X));
@constraint(m, tr(A1 * X) == b1);
@constraint(m, tr(A2 * X) == b2);
@constraint(m, psdcon, X in JuMP.PSDCone()); # this is used to enforce a PsdConeSquare constraint
status = JuMP.optimize!(m);

@testset "SDP Square - JuMP" begin
    @test isapprox(JuMP.objective_value(m), 13.9022, atol = 1e-3)
    @test JuMP.termination_status(m) == MathOptInterface.Success
    @test JuMP.primal_status(m) == MathOptInterface.FeasiblePoint
end

# ---------------------------
# Semidefinite Program - Triangle
# ---------------------------
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

    m = Model(with_optimizer(COSMO.Optimizer));
    @variable(m, X[1:3,1:3], PSD);
    @objective(m, Min, tr(C * X));
    @constraint(m, tr(A1 * X) == b1);
    @constraint(m, tr(A2 * X) == b2);
    status = JuMP.optimize!(m);

@testset "SDP Triangle - JuMP" begin
    @test isapprox(JuMP.objective_value(m), 13.9022, atol = 1e-3)
    @test JuMP.termination_status(m) == MathOptInterface.Success
    @test JuMP.primal_status(m) == MathOptInterface.FeasiblePoint
end

nothing