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

function admmStep(x, s, μ, ν, x_, s_, ls, k, F, q, b, K, ρ, α, σ, m, n)
  @. ls[1:n] = σ.*x-q
  @. ls[n+1:end] = b-s+μ./ρ
  # k = cg(F, RHS - F*[x_; ν], 10) + [x_; ν]
  ls = F \ ls
  #deconstruct solution vector k = [x(n+1);ν(n+1)]
  @. x_ = k[1:n]
  @. ν = k[n+1:end]
  # Projection steps
  @. x = α*x_ + (1.0-α)*x
  @. s_ = s - (ν+μ)./ρ
  @. s_ = α*s_ + (1.0-α)*s
  @. s = s_ + μ./ρ
  Projections.projectCompositeCone!(s, K)

  # update dual variable μ
  @. μ = μ + ρ.*(s_ - s)
  nothing
end
# SOLVER ROUTINE
# -------------------------------------
  function solve(P,q,A,b,K::OSSDPTypes.Cone,settings::OSSDPTypes.OSSDPSettings)

    # create workspace variables
    ws = WorkSpace(Problem(P,q,A,b,K),ScaleMatrices())

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

    tic()
    startTime = time()
    xtilde = zeros(ws.p.n)
    stilde = zeros(ws.p.m)

    # MAIN ADMM LOOP
    # F = copy(ws.p.F)
    x_ = similar(ws.x)
    s_ = similar(ws.s)
    const n = ws.p.n
    const m = ws.p.m
    ls = zeros(n + m)
    k = similar()
    for iter = 1:settings.max_iter
      admmStep(
        ws.x, ws.s, ws.μ, ws.ν,
        x_, s_, ls, k,
        ws.p.F, ws.p.q, ws.p.b, K, ws.p.ρVec,
        settings.alpha, settings.sigma,
        m, n
      )

      # compute residuals (based on optimality conditions of the problem) to check for termination condition
      # compute them every {settings.checkTermination} step
      mod(iter,settings.checkTermination)  == 0 && ((r_prim,r_dual) = calculateResiduals(ws,settings))

      # compute deltas
      # δx = xNew - xPrev
      # δy = yNew - yPrev

      # print iteration steps
      if settings.verbose && mod(iter,settings.checkTermination)  == 0
        # update cost
        cost = ws.sm.cinv*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]
        printIteration(settings,iter,cost,r_prim,r_dual)
      end


      # check convergence with residuals every {settings.checkIteration} step
      if mod(iter,settings.checkTermination) == 0
        if hasConverged(ws,settings,r_prim,r_dual)
          status = :solved
          break
        end
      end

      # adapt rhoVec if enabled
      if settings.adaptive_rho && (mod(iter,settings.adaptive_rho_interval) == 0) && (settings.adaptive_rho_interval > 0)
        adaptRhoVec!(ws,settings)
      end

      if settings.timelimit !=0 &&  (time() - startTime) > settings.timelimit
        status = :TimeLimit
        break
      end

    end #END-ADMM-MAIN-LOOP

    rt = toq()

    # calculate primal and dual residuals
    if iter == settings.max_iter
      r_prim,r_dual = calculateResiduals(ws,settings)
      status = :UserLimit
    end

    if settings.scaling != 0
      reverseScaling!(ws)
      # FIXME: Another cost calculation is not necessary since cost value is not affected by scaling
      cost =  (1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1] #sm.cinv * not necessary anymore since reverseScaling
    end

    # print solution to screen
    settings.verbose && printResult(status,iter,cost,rt)


    # create result object
    result = OSSDPResult(ws.x,ws.s,ws.ν,ws.μ,cost,iter,status,rt,r_prim,r_dual);

    return result,ws;

  end

end




