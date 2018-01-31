module OSSDPTypes

export sdpResult, sdpDebug, problem, sdpSettings,scaleMatrices
# -------------------------------------
# TYPE DEFINITIONS
# -------------------------------------
  type sdpResult
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

  type sdpDebug
    x::Array{Float64,2}
    s::Array{Float64,2}
    λ::Array{Float64,2}
    ν::Array{Float64,2}
    cost::Array{Float64}
  end

  type problem
    P::Array{Float64,2}
    q::Array{Float64,1}
    A::Array{Float64,2}
    b::Array{Float64,1}
    m::Int64
    n::Int64
  end

  type scaleMatrices
    D::Array{Float64,2}
    E::Array{Float64,2}
    scaleMatrices() = new()
  end


  type sdpSettings
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
    # Scaling
    scaling::Int64
    MIN_SCALING::Float64
    MAX_SCALING::Float64


    #constructor
    function sdpSettings(;
      rho=1.0,
      sigma=10e-6,
      alpha=1.6,
      eps_abs=1e-5,
      eps_rel=1e-5,
      eps_prim_inf=1e-4,
      eps_dual_inf=1e-4,
      max_iter=2500,
      verbose=false,
      checkTermination=1,
      scaling=10,
      MIN_SCALING = 1e-04,
      MAX_SCALING = 1e+04
      )
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose,checkTermination,scaling,MIN_SCALING,MAX_SCALING)
    end
  end

  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::sdpResult)
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\n\nx = $(round.(obj.x,3))\ns = $(round.(obj.s,3))\nz = $(round.(obj.z,3))\nμ = $(round.(obj.μ,3))\nλ = $(round.(obj.λ,3))" )
  end
end
