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
    status::String
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
    P::Array{Float64,2}
    q::Array{Float64,1}
    A::Array{Float64,2}
    b::Array{Float64,1}
    m::Int64
    n::Int64
  end

  mutable struct scaleMatrices
    D::Array{Float64,2}
    Dinv::Array{Float64,2}
    E::Array{Float64,2}
    Einv::Array{Float64,2}
    c::Float64
    cinv::Float64
    scaleMatrices() = new(zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),0.,0.)
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
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\n\nx = $(round.(obj.x,3))\ns = $(round.(obj.s,3))\nz = $(round.(obj.z,3))\nμ = $(round.(obj.μ,3))\nλ = $(round.(obj.λ,3))" )
  end
end
