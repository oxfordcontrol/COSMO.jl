
# -------------------------------------
# struct DEFINITIONS
# -------------------------------------
  mutable struct Result
    x::Array{Float64}
    s::Array{Float64}
    ν::Array{Float64}
    μ::Array{Float64}
    cost::Float64
    iter::Int64
    status::Symbol
    solverTime::Float64
    setupTime::Float64
    iterTime::Float64
    rPrim::Float64
    rDual::Float64

    function Result()
      return new(Float64[],Float64[],Float64[],Float64[],0.,0,:Unsolved,0.,0.,0.,0.,0.)
    end

    function Result(x,s,ν,μ,cost,iter,status,solverTime,setupTime,iterTime,rPrim,rDual)
      return new(x,s,ν,μ,cost,iter,status,solverTime,setupTime,iterTime,rPrim,rDual)
    end

  end


  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::Result)
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\nSetup Time: $(round.(obj.setupTime*1000,2))ms\nAvg Iter Time: $(round.((obj.iterTime/obj.iter)*1000,2))ms" )
  end


  # product of cones dimensions, similar to SeDuMi
  mutable struct Cone
    # number of zero  components
    f::Int64
    # number of nonnegative components
    l::Int64
    # dimensions of lorentz constraints (if multiple cones than it's an array)
    q::Array{Int64}
    # dimension of positive semidefinite (psd) constraints
    s::Array{Int64}

    #constructor

    function Cone()
      return new(0,0,Int64[],Int64[])
    end

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
    K::Cone
    ρVec::Array{Float64,1}
    Info::Info
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
      p::Problem
      sm::ScaleMatrices
      x::Vector{Float64}
      s::Vector{Float64}
      ν::Vector{Float64}
      μ::Vector{Float64}
      #constructor
    function WorkSpace(p::Problem,sm::ScaleMatrices)
      m = p.m
      n = p.n
      new(p,sm,zeros(n),zeros(m),zeros(m),zeros(m))
    end
  end


mutable struct Model
    P::SparseMatrixCSC{Float64,Int64}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int64}
    b::Vector{Float64}
    K::Cone

    function Model()
        return new(spzeros(Float64,1,1), Float64[], spzeros(Float64,1,1),Float64[],Cone())
    end

    # function Model(P::SparseMatrixCSC{Float64,Int64},q::Vector{Float64},A::SparseMatrixCSC{Float64,Int64},b::Vector{Float64},K::Cone)
    #     return new(P,q,A,b,K)
    # end

end

  mutable struct Settings
    rho::Float64
    sigma::Float64
    alpha::Float64
    eps_abs::Float64
    eps_rel::Float64
    eps_prim_inf::Float64
    eps_dual_inf::Float64
    max_iter::Int64
    verbose::Bool
    check_termination::Int64
    check_infeasibility::Int64
    scaling::Int64
    MIN_SCALING::Float64
    MAX_SCALING::Float64
    adaptive_rho::Bool
    adaptive_rho_interval::Int64
    adaptive_rho_tolerance::Float64
    RHO_MIN::Float64
    RHO_MAX::Float64
    RHO_TOL::Float64
    timelimit::Int64
    obj_true::Float64
    obj_true_tol::Float64
    #constructor
    function Settings(;
      rho=0.1,
      sigma=1e-6,
      alpha=1.6,
      eps_abs=1e-4,
      eps_rel=1e-4,
      eps_prim_inf=1e-6,
      eps_dual_inf=1e-4,
      max_iter=2500,
      verbose=false,
      check_termination=40,
      check_infeasibility=40,
      scaling=10,
      MIN_SCALING = 1e-4,
      MAX_SCALING = 1e4,
      adaptive_rho = true,
      adaptive_rho_interval = 40,
      adaptive_rho_tolerance = 5,
      RHO_MIN = 1e-6,
      RHO_MAX = 1e6,
      RHO_TOL = 1e-4,
      timelimit = 0,
      obj_true = NaN,
      obj_true_tol = 1e-3
      )
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose,
          check_termination,check_infeasibility,scaling,MIN_SCALING,MAX_SCALING,adaptive_rho,
          adaptive_rho_interval,adaptive_rho_tolerance,RHO_MIN,RHO_MAX,RHO_TOL,timelimit,obj_true,obj_true_tol)
    end
  end



abstract type AbstractConvexSet end


struct ZeroCone <:AbstractConvexSet
  A::AbstractMatrix{<:Real}
  b::AbstractVector{<:Real}
  dim::Int
  project!::Function
end

struct Box{T <: Real} <: AbstractConvexSet
    A::AbstractMatrix{<:Real}
    b::AbstractVector{<:Real}
    dim::Int
    project!::Function
    l::T
    u::T
end

struct NonNegativeOrthant <:AbstractConvexSet
  A::AbstractMatrix{<:Real}
  b::AbstractVector{<:Real}
  dim::Int
  project!::Function
end

struct SecondOrderCone <:AbstractConvexSet
  A::AbstractMatrix{<:Real}
  b::AbstractVector{<:Real}
  dim::Int
  project!::Function
end

struct PositiveSemidefiniteCone <:AbstractConvexSet
  A::AbstractMatrix{<:Real}
  b::AbstractVector{<:Real}
  dim::Int
  project!::Function
end


