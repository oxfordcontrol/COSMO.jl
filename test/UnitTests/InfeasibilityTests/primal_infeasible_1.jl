

# Primal infeasible test problems
# x >= 0, elements in A >=0, elements in b <0 and s in R+

nn = 1
rng = Random.MersenneTwister(1313)

@testset "Primal infeasible QP problems - Testset 1" begin
  for iii = 1:nn

    # choose size of problem
    n = rand(rng,5:50)
    m = 2*n
    A = sprand(rng,m,n,0.8)*50
    b = -rand(rng,m)*50
    A = [A; -sparse(1.0I,n,n)]
    b = [b;zeros(n)]

    # create dual feasibile problem
    P = Helper.generatePosDefMatrix(n,rng)
    ytrue = rand(rng,m+n,1)*50
    xtrue = rand(rng,n,1)*50
    q = (-P*xtrue -  A'*ytrue)[:]

    constraint = QOCS.Constraint(-A,b,QOCS.Nonnegatives())
    ra = 0.

    settings = QOCS.Settings(max_iter=10000,eps_abs = 1e-5,eps_rel=1e-5)
    model = QOCS.Model()
    assemble!(model,P,q,[constraint])
    res1 = QOCS.optimize!(model,settings);

     @test res1.status == :Primal_infeasible
  end
end
nothing
  # # solve accurately once with mosek s
  # model = Model(solver=MosekSolver())
  # @variable(model, x[1:n])
  # @variable(model, s[1:m+n])

  # @objective(model, Min, 0.5*x'*P*x+q'*x)
  # @constraint(model,A*x+s.==b)
  # @constraint(model,s.>=0)
  # status = JuMP.solve(model)


  # end


