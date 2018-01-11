module OSSDP

using Formatting, Projections
export solveSDP, sdpResult, sdpDebug, sdpSettings

# -------------------------------------
# TYPE DEFINITIONS
# -------------------------------------
  type sdpResult
    x::Array{Float64}
    s::Array{Float64}
    z::Array{Float64}
    λ::Array{Float64}
    μ::Array{Float64}
    cost::Float64
    iter::Int64
    status::String
    solverTime::Float64
  end

  type sdpDebug
    x::Array{Float64,2}
    s::Array{Float64,2}
    z::Array{Float64,2}
    λ::Array{Float64,2}
    μ::Array{Float64,2}
    ν::Array{Float64,2}
    cost::Array{Float64}
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

    #constructor
    function sdpSettings(;rho=1.0,sigma=10e-6,alpha=1.6,eps_abs=1e-3,eps_rel=1e-3,eps_prim_inf=1e-4,eps_dual_inf=1e-4,max_iter=2500,verbose=false)
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose)
    end
  end

  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::sdpResult)
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\n\nx = $(round.(obj.x,3))\ns = $(round.(obj.s,3))\nz = $(round.(obj.z,3))\nμ = $(round.(obj.μ,3))\nλ = $(round.(obj.λ,3))" )
  end




  function isPrimalInfeasible(δy::Array{Float64},A,l::Array{Float64},u::Array{Float64},ϵ_prim_inf::Float64)
    norm_δy = norm(δy,Inf)
    if norm_δy > ϵ_prim_inf
      δy = δy/norm_δy
      # second condition
      if (u'*max.(δy,0) + l'*min.(δy,0) )[1] <= - ϵ_prim_inf*norm_δy
        # first condition
          if norm(A'*δy,Inf) <= ϵ_prim_inf*norm_δy
           return true
          end
        end
      end
      return false
  end


  function isDualInfeasible(δx::Array{Float64},P,A,q::Array{Float64},l::Array{Float64},u::Array{Float64},ϵ_dual_inf::Float64)
    norm_δx = norm(δx,Inf)
    m = size(A,1)
    if norm_δx > ϵ_dual_inf
      if (q'*δx)[1] < - ϵ_dual_inf*norm_δx
        if norm(P*δx,Inf) < ϵ_dual_inf*norm_δx
          Aδx = A * δx
          for i = 1:m
            if ( (u[i] < 1e18) && (Aδx[i] > ϵ_dual_inf*norm_δx) ) || ( (l[i] > -1e18) && (Aδx[i] < - ϵ_dual_inf*norm_δx) )
              return false
            end
          end
          return true
        end
      end
    end
    return false
  end

# SOLVER ROUTINE
# -------------------------------------
  function solveSDP(P::Array{Float64,2},q::Array{Float64,1},A::Array{Float64,2},b::Array{Float64,1},settings::sdpSettings)

    #Load algorithm settings
    σ = settings.sigma
    α = settings.alpha
    ρ = settings.rho
    ϵ_abs = settings.eps_abs
    ϵ_rel = settings.eps_rel
    ϵ_prim_inf = settings.eps_prim_inf
    ϵ_dual_inf = settings.eps_dual_inf

    # instantiate variables
    iter = 0
    status = "unsolved"
    r_prim = 100
    r_dual = 100

    # determine size of decision variables
    # n: r^2 since we are working with vectorized matrixes of size r
    n = size(q,1)
    r = Int(sqrt(n))
    m = size(b,1)
    xPrev = zeros(n)
    sPrev = zeros(n)
    zPrev = zeros(m)
    λPrev = zeros(n)
    μPrev = zeros(m)
    xNew = zeros(n)
    sNew = zeros(n)
    zNew = zeros(m)
    λNew = zeros(n)
    μNew = zeros(m)
    ν = zeros(m)
    cost = Inf

    # create debugging variables
    xArr = zeros(settings.max_iter,n)
    sArr = zeros(settings.max_iter,n)
    zArr = zeros(settings.max_iter,m)
    λArr = zeros(settings.max_iter,n)
    μArr = zeros(settings.max_iter,m)
    costArr = zeros(settings.max_iter)
    νArr = zeros(settings.max_iter,m)

    # TODO: Print information still true?
    # print information about settings to the screen
    println("-"^50 * "\n" * " "^8 * "ADMM-SDP Solver in pure Julia\n" * " "^18 * "Michael Garstka\n"  * " "^8 * "University of Oxford, January 2018\n" * "-"^50 * "\n")
    println("Problem: variable vec(X) size: n = $(n), constraints m = $(m)")
    println("Settings: ϵ_abs = $(ϵ_abs), ϵ_rel = $(ϵ_rel),\n" * " "^10 * "ϵ_prim_inf = $(ϵ_prim_inf), ϵ_dual_inf = $(ϵ_dual_inf),\n" * " "^10 * "ρ = $(ρ), σ = $(σ), α = $(α),\n" * " "^10 * "max_iter = $(settings.max_iter)\n\n")

    tic()

    # KKT matrix M
    # FIXME: Correct representation of P
    P=0.0
    M = [P+σ*eye(n) A';A -(1/ρ)*eye(m)]
    M = sparse(M)

    # Do LDLT Factorization: A = LDL^T
    F = ldltfact(M)

    for iter = 1:1:settings.max_iter

      # assign previous variables
      xPrev = xNew
      sPrev = sNew
      zPrev = zNew

      # construct right hand side
      RHS = [-q+σ*sPrev-λPrev; b-(1/ρ)*μPrev]

      #solve linear system M*k = b with help of factorization matrix
      k = F\RHS

      #deconstruct solution vector k = [xt_(k+1);ν_(k+1)]
      xt = k[1:n]
      ν = k[n+1:end]
      zt = b + (1/ρ) * (ν - μPrev)

      # Projection steps and relaxation
      # TODO: Find out why and where relaxation with α makes sense
      xNew = α*xt + (1-α)*xPrev

      #TODO: SCS uses approximate projection (see Paper)
      sNew = Projections.sdcone( xNew + (1/σ)*λPrev,r)
      #zNew = α*zt + (1-α)*zPrev
      zNew = zt

      # update dual variables
      λNew = λPrev + σ*(xNew - sNew)
      μNew = μPrev + ρ*(zNew - b)

      # update cost
      cost = (1/2 * xNew'*P*xNew + q'*xNew)[1]

      # compute residuals to check for termination condition
      # TODO: Correct residuals?
      r_prim = norm(A*xNew - zNew,Inf)
      #r_dual = norm(P*xNew + q + A'*μNew + λNew,Inf)
      r_dual = norm(σ*(sNew - sPrev),Inf)

      # store variables
      xArr[iter,:] = xNew
      sArr[iter,:] = sNew
      zArr[iter,:] = zNew
      λArr[iter,:] = λNew
      μArr[iter,:] = μNew
      costArr[iter] = cost
      νArr[iter,:] = ν
      # compute deltas
      # δx = xNew - xPrev
      # δy = yNew - yPrev

      # print iteration steps
      if settings.verbose
        if iter == 1
          println("Iter:\tObjective:\tPrimal Res:\tDual Res:")
        end
        if mod(iter,100) == 0 || iter == 1 || iter == 2 || iter == settings.max_iter
          printfmt("{1:d}\t{2:.4e}\t{3:.4e}\t{4:.4e}\n", iter,cost,r_prim,r_dual)
       end
      end

      # if isPrimalInfeasible(δy,A,l,u,ϵ_prim_inf)
      #     status = "primal infeasible"
      #     cost = Inf
      #     xNew = NaN*ones(n,1)
      #     yNew = NaN*ones(m,1)
      #     break
      # end

      # if isDualInfeasible(δx,P,A,q,l,u,ϵ_dual_inf)
      #     status = "dual infeasible"
      #     cost = -Inf
      #     xNew = NaN*ones(n,1)
      #     yNew = NaN*ones(m,1)
      #     break
      # end

      # check convergence with residuals
      # TODO: Check convergence condition for SDP case
      ϵ_prim = ϵ_abs + ϵ_rel * max(norm(A*xNew,Inf), norm(zNew,Inf) )
      ϵ_dual = ϵ_abs + ϵ_rel * max(norm(P*xNew,Inf), norm(A'*μNew,Inf), norm(q,Inf) )
      if ( r_prim < ϵ_prim && r_dual < ϵ_dual)
        if settings.verbose
          printfmt("{1:d}\t{2:.4e}\t{3:.4e}\t{4:.4e}\n", iter,cost,r_prim,r_dual)
        end
        status = "solved"
        break
      end
    end

    # print solution to screen
    rt = toq()
    println("\n\n" * "-"^50 * "\nRESULT: Status: $(status)\nTotal Iterations: $(iter)\nOptimal objective: $(round.(cost,4))\nRuntime: $(round.(rt,3))s ($(round.(rt*1000,2))ms)\n" * "-"^50 )

    # create result object
    #TODO: Change result object
    result = sdpResult(xNew,sNew,zNew,λNew,μNew,cost,iter,status,rt);

    dbg = sdpDebug(xArr,sArr,zArr,λArr,μArr,νArr,costArr)

    return result,dbg;

  end

end




