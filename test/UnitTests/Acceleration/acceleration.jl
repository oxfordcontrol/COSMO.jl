using COSMO, Test, Random, LinearAlgebra

rng = Random.MersenneTwister(232311)
ea = COSMO.EmptyAccelerator{Float64}()

# x = rand(rng, 2)
# g = rand(rng, 2)

# # @testset "Acceleration" begin
#   # Check empty accelerator
#   @test COSMO.update_history!(ea, x, g) == nothing
#   @test COSMO.accelerate!(g, x, ea)

#   # Check BLAS.gesv! wrapper
#   A = [1. 2; 0 1]
#   A_f = [1. 2; 0 0]
#   b = rand(rng, 2)
#   @test COSMO._gesv!(A, b) == 1
#   @test COSMO._gesv!(A_f, b) == -1

#   # Check Anderson Accelerator
#   @test_throws DomainError COSMO.AndersonAccelerator{Float64}(0)
#   @test_throws DomainError COSMO.AndersonAccelerator{Float64}(1, 0)

  mem = 10
  dim = 2
  aa = COSMO.AndersonAccelerator{Float64}(dim, mem = mem, is_type1 = true);
  xmin = -10.
  xmax = 10.
  # create some reference iteration data
  # x_iter = xmin .+ randn(rng, dim, mem + 5) * (xmax - xmin)
  # g_iter = xmin .+ randn(rng, dim, mem + 5) * (xmax - xmin)
  x_iter = [[6.; 6.] [-1.8; -4.8] [1.56; 2.04] [-1.2; -1.344] [-0.2736; 0.9888] [-1.30272; -0.0336] [-1.10093; 0.788352] [-1.151338; 0.502886] [-1.50708;  0.807452] [-1.68731; 0.74276] [-1.72058; 0.863831] [-1.80653; 0.859581] [-1.83836; 0.912002]   ]
  g_iter = [[-1.8; -4.8] [1.56; 2.04] [-1.2; -1.344] [-0.2736; 0.9888] [-1.30272; -0.0336] [-1.10093; 0.788352] [-1.151338; 0.502886] [-1.50708;  0.807452] [-1.68731; 0.74276] [-1.72058; 0.863831] [-1.80653; 0.859581] [-1.83836; 0.912002] [-1.88255; 0.920616]  ]


  f_iter = x_iter - g_iter
  Xref = zeros(dim, mem + 5 )
  Gref = zeros(dim, mem + 5 )
  Fref = zeros(dim, mem + 5 )

  # pre-calculate the difference matrices
  for j = 2:size(x_iter, 2)-1
    Xref[:, j] = x_iter[:, j] - x_iter[:, j - 1]
    Gref[:, j] = g_iter[:, j] - g_iter[:, j - 1]
    Fref[:, j] = f_iter[:, j] - f_iter[:, j - 1]
  end
  Xref[:, 1] = x_iter[:, 1]
  Gref[:, 1] = g_iter[:, 1]
  Fref[:, 1] = f_iter[:, 1]
  ref_k = 1


 @testset "Acceleration" begin


  #_gesv
  # @test COSMO._gesv!(rand(rng, 2, 2), rand(rng, 2) ) == 1
  # @test COSMO._gesv!([1. 2; 0 0], rand(rng, 2) ) == -1

 # update the history for some iterations
  for i = 1:3
    global ref_k
    COSMO.update_history!(aa, g_iter[:, ref_k], x_iter[:, ref_k])
    ref_k += 1
  end

  @test aa.iter == 3
  @test aa.x_last == x_iter[:, 3]
  @test aa.g_last == g_iter[:, 3]
  @test aa.f_last == f_iter[:, 3]
  @test aa.X == [Xref[:, 1:3] zeros(dim, mem - 3)]
  @test aa.G == [Gref[:, 1:3] zeros(dim, mem - 3)]
  @test aa.F == [Fref[:, 1:3] zeros(dim, mem - 3)]
  @test norm(aa.M[1:3, 1:3] - Xref[:, 1:3]'  * Fref[:, 1:3], Inf) < 1e-10

  # accelerate iterate
  x_acc = copy(x_iter[:, 3])
  g_acc = copy(g_iter[:, 3])
  eta = (Xref[:, 1:3]' * Fref[:, 1:3]) \ (Xref[:, 1:3]' * f_iter[:, 3])
  g_ref = g_acc - Gref[:, 1:3] * eta

  COSMO.accelerate!(g_acc, x_acc, aa)
  @test g_acc ==  g_ref

  #... update history further
  for i = 4:10
    global ref_k
    COSMO.update_history!(aa, g_iter[:, ref_k], x_iter[:, ref_k])
    ref_k += 1
  end

  @test aa.iter == 10
  @test aa.x_last == x_iter[:, 10]
  @test aa.g_last == g_iter[:, 10]
  @test aa.f_last == f_iter[:, 10]
  @test aa.X == Xref[:, 1:10]
  @test aa.G == Gref[:, 1:10]
  @test aa.F == Fref[:, 1:10]
  @test norm(aa.M - Xref[:, 1:10]'  * Fref[:, 1:10], Inf) < 1e-10

  # ... and update even further exceeding the history memory
  for i = 11:12
    global ref_k
    COSMO.update_history!(aa, g_iter[:, ref_k], x_iter[:, ref_k])
    ref_k += 1
  end
  @test aa.iter == 12
  @test aa.X[:, 1:2] == Xref[:, 11:12]
  @test aa.X[:, 3:10] == Xref[:, 3:10]
  @test aa.G[:, 1:2] == Gref[:, 11:12]
  @test aa.F[:, 1:2] == Fref[:, 11:12]


  # accelerate iterate
  x_acc = copy(x_iter[:, 12])
  g_acc = copy(g_iter[:, 12])
  eta = (Xref[:, 3:12]' * Fref[:, 3:12]) \ (Xref[:, 3:12]' * f_iter[:, 12])
  g_ref = g_acc - Gref[:, 1:3] * eta

  COSMO.accelerate!(g_acc, x_acc, aa)
  @test g_acc ==  g_ref

end





# end
nothing
