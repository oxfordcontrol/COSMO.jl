module OSSDPTypes
export OSSDPResult, Problem, OSSDPSettings, ScaleMatrices, Cone, WorkSpace
# -------------------------------------
# struct DEFINITIONS
# -------------------------------------
  struct OSSDPResult
    x::Array{Float64}
    s::Array{Float64}
    ν::Array{Float64}
    μ::Array{Float64}
    cost::Float64
    iter::Int64
    status::Symbol
    solverTime::Float64
    setupTime::Float64
    avgIterTime::Float64
    rPrim::Float64
    rDual::Float64
  end


  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::OSSDPResult)
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\nSetup Time: $(round.(obj.setupTime*1000,2))ms\nAvg Iter Time: $(round.(obj.avgIterTime*1000,2))ms" )
  end


  # product of cones dimensions, similar to SeDuMi
  struct Cone
    # number of zero  components
    f::Int64
    # number of nonnegative components
    l::Int64
    # dimensions of lorentz constraints (if multiple cones than it's an array)
    q::Array{Int64}
    # dimension of positive semidefinite (psd) constraints
    s::Array{Int64}

    #constructor
    function Cone(f::Int64,l::Int64,q,s)
      (f < 0 || l < 0) && error("Negative values are not allowed.")
      (length(q) == 1 && q[1] == 0) && (q = [])
      (length(s) == 1 && s[1] == 0) && (s = [])
      (length(q) > 0 && minimum(q) <= 0) && error("Cone dimensions in K.q have to be positive integers.")
      (length(s) > 0 && minimum(s) <= 0) && error("Cone dimension in K.s have to be positive integers.")
      new(f,l,q,s)
    end
  end
 mutable struct Info
    rho_updates::Array{Float64,1}
  end

  mutable struct Problem
    P::SparseMatrixCSC{Float64,Int64}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int64}
    b::Vector{Float64}
    m::Int64
    n::Int64
    K::OSSDPTypes.Cone
    ρVec::Array{Float64,1}
    Info::OSSDPTypes.Info
    F
    #constructor
    function Problem(P,q,A,b,K)
      # check dimensions
      m = size(A,1)
      n = size(A,2)
      (size(P,1) != n || size(P,2) != n) && error("Dimensions of P and A dont match.")
      (size(q,1) != n || size(q,2) != 1) && error("Dimensions of P and q dont match.")
      (size(b,1) != m || size(b,2) != 1) && error("Dimensions of A and b dont match.")
      (in(NaN,b) || in(Inf,b) || in(-Inf,b)) && error("b must not contain Inf or NaN.")
      (P != P') && error("P must be symmetric!")
      # Make sure problem data is in sparse format
      typeof(P) != SparseMatrixCSC{Float64,Int64} && (P = sparse(P))
      typeof(A) != SparseMatrixCSC{Float64,Int64} && (A = sparse(A))
      typeof(b) == SparseVector{Float64,Int64} && (b = full(b))
      typeof(q) == SparseVector{Float64,Int64} && (q = full(q))

      # check that number of cone variables provided in K add up
      isempty(K.q) ? nq = 0 :  (nq = sum(K.q) )
      isempty(K.s) ? ns = 0 :  (ns = sum(K.s) )
      (K.f + K.l + nq + ns ) != m && error("Problem dimension doesnt match cone sizes provided in K.")
      # FIXME: prevent copying of data for better performance
      new(copy(P),copy(q),copy(A),copy(b),m,n,K,[0.],Info([0.]),0)
      #new(P,q,A,b,m,n,K) using this seems to change the input data of main solveSDP function
    end
  end

  mutable struct ScaleMatrices
    D::SparseMatrixCSC{Float64,Int64}
    Dinv::SparseMatrixCSC{Float64,Int64}
    E::SparseMatrixCSC{Float64,Int64}
    Einv::SparseMatrixCSC{Float64,Int64}
    c::Float64
    cinv::Float64
    ScaleMatrices() = new(spzeros(1,1),spzeros(1,1),spzeros(1,1),spzeros(1,1),1.,1.)
  end

  mutable struct WorkSpace
      p::OSSDPTypes.Problem
      sm::OSSDPTypes.ScaleMatrices
      x::Vector{Float64}
      s::Vector{Float64}
      ν::Vector{Float64}
      μ::Vector{Float64}
      #constructor
    function WorkSpace(p::OSSDPTypes.Problem,sm::OSSDPTypes.ScaleMatrices)
      m = p.m
      n = p.n
      new(p,sm,zeros(n),zeros(m),zeros(m),zeros(m))
    end
  end

  mutable struct OSSDPSettings
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
    checkInfeasibility::Int64
    scaling::Int64
    MIN_SCALING::Float64
    MAX_SCALING::Float64
    avgFunc::Function
    scaleFunc::Int64
    adaptive_rho::Bool
    adaptive_rho_interval::Int64
    adaptive_rho_tolerance::Float64
    RHO_MIN::Float64
    RHO_MAX::Float64
    RHO_TOL::Float64
    timelimit::Int64
    objTrue::Float64
    objTrueTOL::Float64
    #constructor
    function OSSDPSettings(;
      rho=0.1,
      sigma=1e-6,
      alpha=1.6,
      eps_abs=1e-4,
      eps_rel=1e-4,
      eps_prim_inf=1e-6,
      eps_dual_inf=1e-4,
      max_iter=2500,
      verbose=false,
      checkTermination=1,
      checkInfeasibility=40,
      scaling=10,
      MIN_SCALING = 1e-4,
      MAX_SCALING = 1e4,
      avgFunc = mean,
      scaleFunc = 2,
      adaptive_rho = false,
      adaptive_rho_interval = 40,
      adaptive_rho_tolerance = 5,
      RHO_MIN = 1e-6,
      RHO_MAX = 1e6,
      RHO_TOL = 1e-4,
      timelimit = 0,
      objTrue = NaN,
      objTrueTOL = 1e-3
      )
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose,
          checkTermination,checkInfeasibility,scaling,MIN_SCALING,MAX_SCALING,avgFunc,scaleFunc,adaptive_rho,
          adaptive_rho_interval,adaptive_rho_tolerance,RHO_MIN,RHO_MAX,RHO_TOL,timelimit,objTrue,objTrueTOL)
    end
  end

end
