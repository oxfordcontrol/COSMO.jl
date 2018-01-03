module OSQPSolver
using Formatting
export solveOSQP, qpResult, qpSettings, test

# TODO: use vec() and reshape() to perform the vectorize and matrixify operations
# TODO: rename OSQP specific objects and function names

# -------------------------------------
# HELPER FUNCTIONS
# -------------------------------------
  # compute projection of x onto a box defined by l and u
  function project_box(x::Array{Float64},l::Array{Float64},u::Array{Float64})
    return min.( max.(x,l), u)
  end

   function project_sdcone(x::Array{Float64},n::Int64)
    # recreate original matrix from input vector
    X = reshape(X,n,n)
    X = X./2
    X = X+X'
    # compute eigenvalue decomposition
    F = eigfact(X)
    Λ = diag(F[:values])
    Q = F[:vectors]
    # set negative eigenvalues to 0
    Xp = Q*max(Λ,0)*Q'
    return vec(Xp)
  end

# -------------------------------------
# TYPE DEFINITIONS
# -------------------------------------
  type qpResult
    x::Array{Float64}
    y::Array{Float64}
    cost::Float64
    iter::Int64
    status::String
    solverTime::Float64
  end
  # Redefinition of the show function that fires when the object is called
  function Base.show(io::IO, obj::qpResult)
    println(io,"\nRESULT: \nTotal Iterations: $(obj.iter)\nCost: $(round.(obj.cost,2))\nStatus: $(obj.status)\nSolve Time: $(round.(obj.solverTime*1000,2))ms\n\nx = $(round.(obj.x,3))\ny = $(round.(obj.y,3))" )
  end


  type qpSettings
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
    function qpSettings(;rho=0.1,sigma=10e-6,alpha=1.5,eps_abs=1e-3,eps_rel=1e-3,eps_prim_inf=1e-5,eps_dual_inf=1e-5,max_iter=2500,verbose=false)
        new(rho,sigma,alpha,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,max_iter,verbose)
    end
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
  function solveOSQP(P,q::Array{Float64},A,l::Array{Float64},u::Array{Float64},b::Array{Float64},settings::qpSettings)

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
    m = size(A,1)
    n = size(A,2)
    x = zeros(n)
    y = zeros(m)
    z = zeros(m)

    xhat = copy(x)
    zhat = copy(z)
    cost = Inf

    # create initial guess vectors
    xPrev = zeros(n,1)
    yPrev = zeros(m,1)
    zPrev = zeros(m,1)
    nuNew = zeros(m,1)
    xNew = zeros(n,1)
    yNew = zeros(m,1)
    zNew = zeros(m,1)

    # print information about settings to the screen
    println("-"^50 * "\n" * " "^8 * "ADMM-SDP Solver in pure Julia\n" * " "^18 * "Michael Garstka\n"  * " "^8 * "University of Oxford, October 2017\n" * "-"^50 * "\n")
    println("Problem:  variables n = $(n), constraints m = $(m)")
    println("Settings: ϵ_abs = $(ϵ_abs), ϵ_rel = $(ϵ_rel),\n" * " "^10 * "ϵ_prim_inf = $(ϵ_prim_inf), ϵ_dual_inf = $(ϵ_dual_inf),\n" * " "^10 * "ρ = $(ρ), σ = $(σ), α = $(α),\n" * " "^10 * "max_iter = $(settings.max_iter)\n\n")

    tic()

    # KKT matrix M
    M = [P+σ*eye(n) A';A -1/ρ*eye(m)]
    M = sparse(M)

    # Do LDLT Factorization: A = LDL^T
    F = ldltfact(M)

    for iter = 1:1:settings.max_iter

      # assign previous variables
      xPrev = xNew
      yPrev = yNew
      zPrev = zNew

      # construct right hand side
      RHS = [σ*xPrev-q; b-zPrev + 1/ρ*yPrev]

      #solve linear system M*k = b with help of factorization matrix
      k = F\RHS

      #deconstruct solution vector k = [xt_(k+1);nu_(k+1)]
      xt = k[1:n]
      nuNew = k[n+1:end]
      zt = zPrev + 1/ρ * (nuNew - yPrev)

      # Projection steps
      xNew = α*xt + (1-α)*xPrev
      zNew = project_sdcone( α*zt + (1-α)*zPrev + 1/ρ*yPrev,n)

      # update dual variable
      yNew = yPrev + ρ* (α*zt + (1-α)*zPrev - zNew)

      # update cost
      cost = (1/2 * xNew'*P*xNew + q'*xNew)[1]

      # compute residuals to check for termination condition
      r_prim = norm(A*xNew - zNew,Inf)
      r_dual = norm(P*xNew + q + A'*yNew,Inf)

      # compute deltas
      δx = xNew - xPrev
      δy = yNew - yPrev

      # print iterations steps
      # TODO: Might be slow
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
      ϵ_prim = ϵ_abs + ϵ_rel * max(norm(A*xNew,Inf), norm(zNew,Inf) )
      ϵ_dual = ϵ_abs + ϵ_rel * max(norm(P*xNew,Inf), norm(A'*yNew,Inf), norm(q,Inf) )
      if ( r_prim < ϵ_prim && r_dual < ϵ_dual  )
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
    result = qpResult(xNew,yNew,cost,iter,status,rt);

    return result;

  end

end




