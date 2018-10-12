
# -------------------------------------
# struct DEFINITIONS
# -------------------------------------

  mutable struct ResultTimes
    solverTime::Float64
    setupTime::Float64
    graphTime::Float64
    factorTime::Float64
    iterTime::Float64
    projTime::Float64
    postTime::Float64

    function ResultTimes()
      return new(0.,0.,0.,0.,0.,0.,0.)
    end

    function ResultTimes(sT::Float64,seT::Float64,gT::Float64,fT::Float64,iT::Float64,pT::Float64,poT::Float64)
      return new(sT,seT,gT,fT,iT,pT,poT)
    end

  end

  struct ResultInfo
    rPrim::Float64
    rDual::Float64
  end
  """
      Result

  Object returned by the QOCS solver after calling `optimize!(model,settings)`. It has the following fields:


  Fieldname | Type | Description
  ---  | --- | ---
  x | Vector{Float64}| Primal variable
  y | Vector{Float64}| Dual variable
  s | Vector{Float64}| (Primal) set variable
  objVal | Float64 | Objective value
  iter | Int64 | Number of iterations
  status | Symbol | Solution status
  info | QOCS.ResultInfo | Struct with more information
  times | QOCS.ResultTimes | Struct with several measured times
  """
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

  function Base.show(io::IO, obj::Result)
    print(io,">>> QOCS - Results\nStatus: $(obj.status)\nIterations: $(obj.iter)\nOptimal Objective: $(round.(obj.objVal,digits=2))\nRuntime: $(round.(obj.times.solverTime*1000,digits=2))ms\nSetup Time: $(round.(obj.times.setupTime*1000,digits=2))ms\n")
    obj.times.iterTime != NaN && print("Avg Iter Time: $(round.((obj.times.iterTime/obj.iter)*1000,digits=2))ms")
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

  struct ScaleMatrices{Tf}
    D::Union{ UniformScaling{Bool}, Diagonal{Tf,Array{Tf,1}}}
    Dinv::Union{UniformScaling{Bool},Diagonal{Tf,Array{Tf,1}} }
    E::Union{UniformScaling{Bool},Diagonal{Tf,Array{Tf,1}} }
    Einv::Union{UniformScaling{Bool},Diagonal{Tf,Array{Tf,1}} }
    c::Base.RefValue{Tf}
    cinv::Base.RefValue{Tf}
  end

  ScaleMatrices()    = ScaleMatrices(Float64)
  ScaleMatrices(m,n) = ScaleMatrices(Float64,m,n)

  function ScaleMatrices(T::Type)
      ScaleMatrices(I,I,I,I,Base.RefValue{T}(one(T)),Base.RefValue{T}(one(T)))
  end

  function ScaleMatrices(T::Type,m,n)
    @assert T <: AbstractFloat
    D    = Diagonal(ones(T,n))
    Dinv = Diagonal(ones(T,n))
    E    = Diagonal(ones(T,m))
    Einv = Diagonal(ones(T,m))
    c    = Base.RefValue{T}(one(T))
    cinv = Base.RefValue{T}(one(T))
    ScaleMatrices(D,Dinv,E,Einv,c,cinv)
  end

mutable struct Flags
  FACTOR_LHS::Bool
  INFEASIBILITY_CHECKS::Bool
  REVERSE_SCALE_PROBLEM_DATA::Bool
  Flags() = new(true,true,true)

end



"""
    Model()

Initializes an empty QOCS model that can be filled with problem data using `assemble!(model,P,q,constraints)`.
"""
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
    M::SparseMatrixCSC{Float64,Int64}
    flags::QOCS.Flags

    function Model()
        return new(spzeros(Float64,1,1), Float64[], spzeros(Float64,1,1),Float64[],AbstractConvexSet[],Cone(),Float64[],Float64[],0,0,0,spzeros(Float64,1,1),Flags())
    end
end

  mutable struct Workspace
      p::QOCS.Model
      sm::ScaleMatrices
      x::Vector{Float64}
      s::Vector{Float64}
      ν::Vector{Float64}
      μ::Vector{Float64}
      ρ::Float64
      ρVec::Array{Float64,1}
      Info::Info
      times::ResultTimes
      #constructor
    function Workspace(p::QOCS.Model,sm::ScaleMatrices)
      m = p.m
      n = p.n
      ws = new(p,sm,zeros(n),zeros(m),zeros(m),zeros(m),0.,Float64[],Info([0.]),ResultTimes())
      # hand over warm starting variables
      length(p.x0) == n && (ws.x = p.x0)
      length(p.y0) == m && (ws.μ = -p.y0)
      return ws
    end
  end

  """
      Settings(;arg=val)

  Creates a QOCS settings object that is used to pass user settings to the solver.

  Argument | Description | Values (default)
  --- | --- | ---
  rho | ADMM rho step | 0.1
  sigma | ADMM sigma step | 1e-6.
  alpha | Relaxation parameter | 1.6
  eps_abs | Absolute residual tolerance | 1e-4
  eps_rel | Relative residual tolerance | 1e-4
  eps_prim_inf | Primal infeasibility tolerance | 1e-4
  eps_dual_inf | Dual infeasibility tolerance | 1e-4
  max_iter | Maximum number of iterations | 2500
  verbose | Verbose printing | false
  verbose_timing | Verbose timing | false
  check_termination | Check termination interval | 40
  check_infeasibility | Check infeasibility interval | 40
  scaling | Number of scaling iterations | 10
  adaptive_rho | Automatic adaptation of step size parameter | true
  time_limit | set solver time limit in s | 0
  """
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
    verbose_timing::Bool
    RHO_MIN::Float64
    RHO_MAX::Float64
    RHO_TOL::Float64
    time_limit::Int64
    obj_true::Float64
    obj_true_tol::Float64
    use_lanczos::Bool
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
      verbose_timing = false,
      RHO_MIN = 1e-6,
      RHO_MAX = 1e6,
      RHO_TOL = 1e-4,
      time_limit = 0,
      obj_true = NaN,
      obj_true_tol = 1e-3,
      use_lanczos = false
      )
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose,
          check_termination,check_infeasibility,scaling,MIN_SCALING,MAX_SCALING,adaptive_rho,
          adaptive_rho_interval,adaptive_rho_tolerance,verbose_timing,RHO_MIN,RHO_MAX,RHO_TOL,time_limit,obj_true,obj_true_tol,
          use_lanczos)
    end
  end

"""
    Constraint(A,b,convexSet,dim=0,indices=0:0)

Creates a QOCS constraint: `Ax + b ∈ convexSet`.

By default the following convex sets are supported: `Zeros`, `Nonnegatives`, `SecondOrderCone`, `PositiveSemidefiniteCone`.

# Examples
```jldoctest
julia> Constraint([1 0;0 1],zeros(2),QOCS.PositiveSemidefiniteCone())
Constraint
Size of A: (2, 2)
ConvexSet: QOCS.PositiveSemidefiniteCone
```

---
The optinal arguments `dim` and `indices` can be used to specify A and b for subparts of variable `x`. If `x` has dimension `dim=4`,
then x[2] and x[3] can be constrained to the zero cone in the following way:


# Examples
```jldoctest
julia> c = Constraint([1 0;0 1],zeros(2),QOCS.Zeros(),4,2:3)
Constraint
Size of A: (2, 4)
ConvexSet: QOCS.Zeros
```
Notice that extra columns of A have been added automatically.
```
julia>Matrix(c.A)
2×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
```
"""
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

    # A = convert(Matrix{Float64},A)
    # b = convert(Vector{Float64},b)

    if indices != 0:0
      (indices.start < 1 || indices.stop < indices.start) && error("The index range for x has to be increasing and nonnegative.")
      dim < indices.stop && error("The dimension of x: $(dim) must be equal or higher than the the stop value of indices: $(indices.stop).")
      Ac = spzeros(size(A,1),dim)
      bc = zeros(dim)
      Ac[:,indices] = A
      bc[indices] = b
      A = Ac
      b = bc
    end
    convexSet.dim = size(A,1)
    #TODO: Change this. I don't think the dimensions of the convex set should be set here, but in their constructor.
    if isa(convexSet, QOCS.PositiveSemidefiniteCone)
      n = Int(sqrt(convexSet.dim))
      convexSet.subspace = zeros(n, n)
    end
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

 function Base.show(io::IO, obj::QOCS.Constraint)
    print(io,"Constraint\nSize of A: $(size(obj.A))\nConvexSet: $(typeof(obj.convexSet))")
  end
