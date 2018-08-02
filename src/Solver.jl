include("./Helper.jl")
include("./Types.jl")
include("./KKT.jl")
include("./Scaling.jl")
include("./Projections.jl")
include("./Residuals.jl")
include("./Parameters.jl")
include("./Infeasibility.jl")
include("./Printing.jl")
include("./Setup.jl")

module OSSDP

using Projections, Scaling, OSSDPTypes, Parameters, Infeasibility, Residuals, Printing, Setup
export solve, OSSDPSettings, Cone #from the Types module


function admmStep!(x, s, μ, ν, x_tl, s_tl, ls, sol, F, q, b, K, ρ, α, σ, m, n)
  # Create right hand side for linear system
  for i=1:n
    ls[i] = σ*x[i]-q[i]
  end
  for i=1:m
    ls[n+i] = b[i]-s[i]+μ[i]/ρ[i]
  end
    sol = F \ ls
  # deconstruct solution vector ls = [x_tl(n+1);ν(n+1)]
  @. x_tl = sol[1:n]
  @. ν = sol[n+1:end]
  # Over relaxation
  @. x = α*x_tl + (1.0-α)*x
  @. s_tl = s - (ν+μ)./ρ
  @. s_tl = α*s_tl + (1.0-α)*s
  @. s = s_tl + μ./ρ
  # Project onto cone K
  Projections.projectCompositeCone!(s, K)
  # update dual variable μ
  @. μ = μ + ρ.*(s_tl - s)
  nothing
end


# SOLVER ROUTINE
# -------------------------------------
  function solve(P,q,A,b,K::OSSDPTypes.Cone,settings::OSSDPTypes.OSSDPSettings)
    runTime_start = time()

    # create workspace variables
    ws = WorkSpace(Problem(P,q,A,b,K),ScaleMatrices())
    P = q = A = b = nothing

    # perform preprocessing steps (scaling, initial KKT factorization)
    tic()
    setup!(ws,settings)
    setupTime = toq()

    # instantiate variables
    iter = 0
    status = :unsolved
    cost = Inf
    r_prim = Inf
    r_dual = Inf


    # print information about settings to the screen
    settings.verbose && printHeader(ws,settings,setupTime)

    timeLimit_start = time()
    #preallocate arrays
    δx = similar(ws.x)
    δy =  similar(ws.μ)
    x_tl = similar(ws.x) # i.e. xTilde
    s_tl = similar(ws.s) # i.e. sTilde
    const n = ws.p.n
    const m = ws.p.m
    ls = zeros(n + m)
    sol = zeros(n + m)

    iter_start = time()

    for iter = 1:settings.max_iter

      @. δx = ws.x
      @. δy = ws.μ
      admmStep!(
        ws.x, ws.s, ws.μ, ws.ν,
        x_tl, s_tl, ls,sol,
        ws.p.F, ws.p.q, ws.p.b, K, ws.p.ρVec,
        settings.alpha, settings.sigma,
        m, n
      )

      # compute deltas for infeasibility detection
      @. δx = ws.x - δx
      @. δy = -ws.μ + δy

      # compute residuals (based on optimality conditions of the problem) to check for termination condition
      # compute them every {settings.checkTermination} step
      mod(iter,settings.checkTermination)  == 0 && ((r_prim,r_dual) = calculateResiduals(ws,settings))


      # check convergence with residuals every {settings.checkIteration} steps
      if mod(iter,settings.checkTermination) == 0
        # update cost
        cost = ws.sm.cinv*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]

        if abs(cost) > 1e20
          status = :Unsolved
          break
        end

        # print iteration steps
        settings.verbose && printIteration(settings,iter,cost,r_prim,r_dual)

        if hasConverged(ws,settings,r_prim,r_dual)
          status = :solved
          break
        end
      end

      # check infeasibility conditions every {settings.checkInfeasibility} steps
      if mod(iter,settings.checkInfeasibility) == 0
        if isPrimalInfeasible(δy,ws,settings)
            status = :primal_infeasible
            cost = Inf
            ws.x .= NaN
            ws.μ .= NaN
            ws.ν .= NaN
            break
        end

        if isDualInfeasible(δx,ws,settings)
            status = :dual_infeasible
            cost = -Inf
            ws.x .= NaN
            ws.μ .= NaN
            ws.ν .= NaN
            break
        end
      end


      # adapt rhoVec if enabled
      if settings.adaptive_rho && (mod(iter,settings.adaptive_rho_interval) == 0) && (settings.adaptive_rho_interval > 0)
        adaptRhoVec!(ws,settings)
      end

      if settings.timelimit !=0 &&  (time() - timeLimit_start) > settings.timelimit
        status = :TimeLimit
        break
      end

    end #END-ADMM-MAIN-LOOP

    iterTime = (time()-iter_start)

    # calculate primal and dual residuals
    if iter == settings.max_iter
      r_prim,r_dual = calculateResiduals(ws,settings)
      status = :UserLimit
    end

    # reverse scaling for scaled feasible cases
    if settings.scaling != 0 && (cost != Inf && cost != -Inf)
      reverseScaling!(ws)
      # FIXME: Another cost calculation is not necessary since cost value is not affected by scaling
      cost =  (1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1] #sm.cinv * not necessary anymore since reverseScaling
    end


    runTime = time() - runTime_start

    # print solution to screen
    settings.verbose && printResult(status,iter,cost,rt)


    # create result object
    result = OSSDPResult(ws.x,ws.s,ws.ν,ws.μ,cost,iter,status,runTime,setupTime,iterTime,r_prim,r_dual);

    return result,ws, δx, -δy;

  end

end




