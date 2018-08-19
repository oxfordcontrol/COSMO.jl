# Unit Test for constraint type

using QOCS, Test, SparseArrays, Random


rng = Random.MersenneTwister(1872381)

A_int = 4
b_int = 2

A_Float = 4.
b_Float = 2.

A_UInt = UInt64(4)
b_UInt = UInt64(2)

A_sparse = sprand(rng,10,2,0.4)
b_sparse = sparse(rand(rng,10))

A_vec = [1 2 3 4]
b_vec = 1

A_mat = rand(rng,10,10)
b_mat = sparse(rand(rng,10,1))

@testset "Constraints" begin

  @testset "Constructors" begin
    @test typeof(QOCS.Constraint(A_int,b_int,QOCS.Zeros())) == QOCS.Constraint
    @test typeof(QOCS.Constraint(A_Float,b_Float,QOCS.Zeros())) == QOCS.Constraint
    @test typeof(QOCS.Constraint(A_UInt,b_UInt,QOCS.Zeros())) == QOCS.Constraint
    @test typeof(QOCS.Constraint(A_sparse,b_sparse,QOCS.Zeros())) == QOCS.Constraint
    @test typeof(QOCS.Constraint(A_vec,b_vec,QOCS.Zeros())) == QOCS.Constraint
    @test typeof(QOCS.Constraint(A_mat,b_mat,QOCS.Zeros())) == QOCS.Constraint
  end

  @testset "Indizes" begin
    A = rand(rng,3,3)
    b = rand(rng,3)
    ind = 3:5
    dim = 10
    cs = QOCS.Constraint(A,b,QOCS.Zeros(),dim,ind)

    @test cs.A[:,ind] == A
    @test cs.b[ind] == b

   end

  @testset "Merge Constraints" begin
    A1 = rand(rng,10,10)
    b1 = rand(rng,10)
    A2 = rand(rng,10,10)
    b2 = rand(rng,10)
    Am = [A1;A2]
    bm = [b1;b2]
    c1 = QOCS.Constraint(A1,b1,QOCS.Zeros())
    c2 = QOCS.Constraint(A2,b2,QOCS.Zeros())
    cArr = [c1;c2]
    cm = QOCS.Constraint(Am,bm,QOCS.Zeros())
    QOCS.mergeConstraints!(cArr)

    @test cArr[1].A == cm.A
    @test cArr[1].b == cm.b
    @test cArr[1].convexSet.dim == cm.convexSet.dim
    @test typeof(cArr[1].convexSet) == typeof(cm.convexSet)


    c1 = QOCS.Constraint(A1,b1,QOCS.Nonnegatives())
    c2 = QOCS.Constraint(A2,b2,QOCS.Nonnegatives())
    cArr = [c1;c2]
    cm = QOCS.Constraint(Am,bm,QOCS.Nonnegatives())
    QOCS.mergeConstraints!(cArr)

    @test cArr[1].A == cm.A
    @test cArr[1].b == cm.b
    @test cArr[1].convexSet.dim == cm.convexSet.dim
    @test typeof(cArr[1].convexSet) == typeof(cm.convexSet)
  end





end
nothing



