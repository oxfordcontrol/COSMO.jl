# unit test for the options factory
using COSMO, Test, LinearAlgebra, SparseArrays, Pardiso

@testset "Options Factory" begin

P = spzeros(1, 1)
q = [1]
c1 = COSMO.Constraint([1.], [-1.], COSMO.Nonnegatives)


iparms = Dict{Int64, Int64}()
iparms[13] = 100
num_threads = 2

settings = COSMO.Settings(kkt_solver = with_options(COSMO.MKLPardisoKKTSolver, iparm = iparms, num_threads = num_threads, msg_level_on = false))
model = COSMO.Model()
assemble!(model, P, q, c1, settings = settings)
model.ρvec = settings.rho * ones(size(model.p.A, 1))
COSMO._make_kkt_solver!(model)

  @test typeof(model.kkt_solver) <: COSMO.MKLPardisoKKTSolver
  # @test Pardiso.get_nprocs(model.kkt_solver.ps) == num_threads
  @test Pardiso.get_iparm(model.kkt_solver.ps, 13) == iparms[13]



  ## check that this works via the MathOptInterface wrapper as well
  optimizer = COSMO.Optimizer(time_limit = 10., kkt_solver = with_options(COSMO.MKLPardisoKKTSolver, iparm = iparms, num_threads = num_threads, msg_level_on = false))
  @test model.settings.kkt_solver == optimizer.inner.settings.kkt_solver

  # assemble model when only KKT solver type is provided
  settings = COSMO.Settings(kkt_solver = COSMO.QdldlKKTSolver)
  model = COSMO.Model()
  assemble!(model, P, q, c1, settings = settings)
  model.ρvec = settings.rho * ones(size(model.p.A, 1))
  COSMO._make_kkt_solver!(model)
  @test typeof(model.kkt_solver) <: COSMO.QdldlKKTSolver

  # try to assemble with wrong type
  struct MockAccelerator end
  @test_throws MethodError COSMO.Settings(kkt_solver = with_options(MockAccelerator, mem = 5))


end


