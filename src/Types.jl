
# -------------------------------------
# struct DEFINITIONS
# -------------------------------------

  struct ResultTimes
    solverTime::Float64
    setupTime::Float64
    iterTime::Float64
  end

  struct ResultInfo
    rPrim::Float64
    rDual::Float64
  end

  struct Result
    x::Array{Float64}
    y::Array{Float64}
    s::Array{Float64}
    objVal::Float64
    iter::Int64
    status::Symbol
    info::ResultInfo
    times::ResultTimes

    # function Result()
    #   return new(Float64[],Float64[],Float64[],0.,0,:Unsolved,ResultInfo(),ResultTimes())
    # end

    # function Result(x,y,s,objVal,iter,status,info,times)
    #   return new(x,y,s,objVal,iter,status,info,times)
    # end

  end

  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::Result)
    print(io,">>> QOCS - Results\nStatus: $(obj.status)\nIterations: $(obj.iter)\nOptimal Objective: $(round.(obj.objVal,digits=2))\nRuntime: $(round.(obj.times.solverTime*1000,digits=2))ms\nSetup Time: $(round.(obj.times.setupTime*1000,digits=2))ms\nAvg Iter Time: $(round.((obj.times.iterTime/obj.iter)*1000,digits=2))ms")
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

  mutable struct ScaleMatrices
    D::SparseMatrixCSC{Float64,Int64}
    Dinv::SparseMatrixCSC{Float64,Int64}
    E::SparseMatrixCSC{Float64,Int64}
    Einv::SparseMatrixCSC{Float64,Int64}
    c::Float64
    cinv::Float64
    ScaleMatrices() = new(spzeros(1,1),spzeros(1,1),spzeros(1,1),spzeros(1,1),1.,1.)
  end

mutable struct Flags
  FACTOR_LHS::Bool
  INFEASIBILITY_CHECKS::Bool
  REVERSE_SCALE_PROBLEM_DATA::Bool
  Flags() = new(true,true,true)

end

mutable struct Model
    P::SparseMatrixCSC{Float64,Int64}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int64}
    b::Vector{Float64}
    convexSets::Array{AbstractConvexSet}
    K::Cone
    x0::Vector{Float64}
    y0::Vector{Float64}
    m::Int64
    n::Int64
    F #LDL Factorization
    flags::QOCS.Flags

    function Model()
        return new(spzeros(Float64,1,1), Float64[], spzeros(Float64,1,1),Float64[],AbstractConvexSet[],Cone(),Float64[],Float64[],0,0,0,Flags())
    end
end

  mutable struct WorkSpace
      p::QOCS.Model
      sm::ScaleMatrices
      x::Vector{Float64}
      s::Vector{Float64}
      ν::Vector{Float64}
      μ::Vector{Float64}
      ρ::Float64
      ρVec::Array{Float64,1}
      Info::Info
      #constructor
    function WorkSpace(p::QOCS.Model,sm::ScaleMatrices)
      m = p.m
      n = p.n
      ws = new(p,sm,zeros(n),zeros(m),zeros(m),zeros(m),0.,Float64[],Info([0.]))
      # hand over warm starting variables
      ws.x = p.x0
      ws.μ = -p.y0
      return ws
    end
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


# TODO: How to handle sparse vs. dense input data
struct Constraint
  A::Union{AbstractMatrix{<:Real},AbstractVector{<:Real}}
  b::AbstractVector{<:Real}
  convexSet::AbstractConvexSet

  function Constraint(A::AbstractMatrix{<:Real},b::AbstractVector{<:Real},convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0)
    size(A,1) != length(b) && error("The dimensions of matrix A and vector b don't match.")
    size(b,2) != 1 && error("Input b must be a vector or a scalar.")
    if dim != 0
      dim < 0 && error("The dimension of x has to be greater than zero.")
    end

    A = convert(Matrix{Float64},A)
    b = convert(Vector{Float64},b)

    if indices != 0:0
      (indices.start < 1 || indices.stop < indices.start) && error("The index range for x has to be increasing and nonnegative.")
      dim < indices.stop && error("The dimension of x: $(dim) must be equal or higher than the the stop value of indices: $(indices.stop).")
      Ac = spzeros(size(A,1),dim)
      bc = zeros(3dim)
      Ac[:,indices] = A
      bc[indices] = b
      A = Ac
      b = bc
    end
    convexSet.dim = size(A,1)
    new(A,b,convexSet)
  end

  function Constraint(A::Real,b::Real,convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0)
    Aarr = zeros(1,1)
    bvec = Vector{Float64}(undef,1)
    Aarr[1] = A
    bvec[1] = b
    Constraint(Aarr,bvec,convexSet,dim,indices)
  end

  Constraint(A::AbstractMatrix{<:Real},b::AbstractMatrix{<:Real},convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0) = Constraint(A,vec(b),convexSet,dim,indices)
  Constraint(A::AbstractMatrix{<:Real},b::Real,convexSet::AbstractConvexSet,dim::Int64=0,indices::UnitRange{Int64}=0:0) = Constraint(A,[b],convexSet,dim,indices)

end



