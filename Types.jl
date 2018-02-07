module OSSDPTypes

export sdpResult, sdpDebug, problem, sdpSettings, scaleMatrices
# -------------------------------------
# struct DEFINITIONS
# -------------------------------------
  struct sdpResult
    x::Array{Float64}
    s::Array{Float64}
    λ::Array{Float64}
    cost::Float64
    iter::Int64
    status::Symbol
    solverTime::Float64
    rPrim::Float64
    rDual::Float64
  end

  struct sdpDebug
    x::Array{Float64,2}
    s::Array{Float64,2}
    λ::Array{Float64,2}
    ν::Array{Float64,2}
    cost::Array{Float64}
  end

  mutable struct problem
    P::SparseMatrixCSC{Float64,Int64}
    q::SparseVector{Float64,Int64}
    A::SparseMatrixCSC{Float64,Int64}
    b::SparseVector{Float64,Int64}
    m::Int64
    n::Int64

    #constructor
    function problem(P,q,A,b)
      # check dimensions
      m = size(A,1)
      n = size(A,2)
      (size(P,1) != n || size(P,2) != n) && error("Dimensions of P and A dont match.")
      (size(q,1) != n || size(q,2) != 1) && error("Dimensions of P and q dont match.")
      (size(b,1) != m || size(b,2) != 1) && error("Dimensions of A and b dont match.")

      # Make sure problem data is in sparse format
      typeof(P) != SparseMatrixCSC{Float64,Int64} && (P = sparse(P))
      typeof(A) != SparseMatrixCSC{Float64,Int64} && (A = sparse(A))
      typeof(b) != SparseVector{Float64,Int64} && (b = sparse(b))
      typeof(q) != SparseVector{Float64,Int64} && (q = sparse(q))
      new(P,q,A,b,m,n)
    end
  end

  mutable struct scaleMatrices
    D::SparseMatrixCSC{Float64,Int64}
    Dinv::SparseMatrixCSC{Float64,Int64}
    E::SparseMatrixCSC{Float64,Int64}
    Einv::SparseMatrixCSC{Float64,Int64}
    c::Float64
    cinv::Float64
    scaleMatrices() = new(spzeros(1,1),spzeros(1,1),spzeros(1,1),spzeros(1,1),1.,1.)
  end


  struct sdpSettings
    rho::Float64
    sigma::Float64
    alpha::Float64
    eps_abs::Float64
    eps_rel::Float64
    eps_prim_inf::Float64
    eps_dual_inf::Float64
    max_iter::Int64
    verbose::Bool
    checkTermination::Int64
    scaling::Int64
    MIN_SCALING::Float64
    MAX_SCALING::Float64


    #constructor
    function sdpSettings(;
      rho=1.0,
      sigma=10.0,
      alpha=1.6,
      eps_abs=1e-6,
      eps_rel=1e-6,
      eps_prim_inf=1e-4,
      eps_dual_inf=1e-4,
      max_iter=2500,
      verbose=false,
      checkTermination=1,
      scaling=10,
      MIN_SCALING = 1e-4,
      MAX_SCALING = 1e4
      )
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose,checkTermination,scaling,MIN_SCALING,MAX_SCALING)
    end
  end

  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::sdpResult)
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\n\nx = $(round.(obj.x,3))\ns = $(round.(obj.s,3))\nν = $(round.(obj.ν,3))\nλ = $(round.(obj.λ,3))" )
  end

  function sparsify!(M)
    M[:,:] = sparse(M)
    return 0
  end
end
