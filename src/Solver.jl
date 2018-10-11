


function admmStep!(x::Vector{Float64}, s::Vector{Float64}, μ::Vector{Float64}, ν::SubArray, x_tl::SubArray, s_tl::Vector{Float64}, ls::Vector{Float64}, sol::Vector{Float64}, F, q::Vector{Float64}, b::Vector{Float64}, ρ::Vector{Float64}, α::Float64, σ::Float64, m::Int64, n::Int64,convexSets::Array{AbstractConvexSet})
  # Create right hand side for linear system
   @. ls[1:n] = σ*x-q
   @. ls[(n+1):end] = b-s+μ/ρ
  sol .= F \ ls
  # deconstruct solution vector ls = [x_tl(n+1);ν(n+1)]
  # x_tl and ν are automatically updated, since they are views on sol

  # Over relaxation
  @. x = α*x_tl + (1.0-α)*x
  @. s_tl = s - (ν+μ)/ρ
  @. s_tl = α*s_tl + (1.0-α)*s
  @. s = s_tl + μ/ρ
  # Project onto cone K
  pTime = @elapsed Projections.projectCompositeCone!(s, convexSets)
  # update dual variable μ
  @. μ = μ + ρ*(s_tl - s)
  return pTime
end



# SOLVER ROUTINE
# -------------------------------------


  """
      optimize!(model,settings)

  Attempts to solve the optimization problem defined in `QOCS.Model` object with the user settings defined in `QOCS.Settings`. Returns a `QOCS.Result` object.
  """
  function optimize!(model::QOCS.Model,settings::QOCS.Settings)
    solverTime_start = time()

    # create scaling variables
    # with scaling    -> uses mutable diagonal scaling matrices
    # without scaling -> uses identity matrices
    sm = (settings.scaling > 0) ? ScaleMatrices(model.m,model.n) : ScaleMatrices()

    # create workspace
    ws = Workspace(model,sm)

    # perform preprocessing steps (scaling, initial KKT factorization)
    ws.times.setupTime = @elapsed setup!(ws,settings);

    # instantiate variables
    numIter = 0
    status = :Unsolved
    cost = Inf
    r_prim = Inf
    r_dual = Inf


    # print information about settings to the screen
    settings.verbose && printHeader(ws,settings)

    timeLimit_start = time()
    #preallocate arrays
    n = ws.p.n
    m = ws.p.m
    δx = zeros(n)
    δy =  zeros(m)
    s_tl = zeros(m) # i.e. sTilde
    n = ws.p.n
    m = ws.p.m
    ls = zeros(n + m)
    sol = zeros(n + m)
    x_tl = view(sol,1:n) # i.e. xTilde
    ν = view(sol,(n+1):(n+m))

    settings.verbose_timing && (iter_start = time())

    for iter = 1:settings.max_iter
      numIter+= 1
      @. δx = ws.x
      @. δy = ws.μ
      ws.times.projTime += admmStep!(
        ws.x, ws.s, ws.μ, ν,
        x_tl, s_tl, ls,sol,
        ws.p.F, ws.p.q, ws.p.b, ws.ρVec,
        settings.alpha, settings.sigma,
        m, n, ws.p.convexSets);

      # compute deltas for infeasibility detection
      @. δx = ws.x - δx
      @. δy = -ws.μ + δy

      # compute residuals (based on optimality conditions of the problem) to check for termination condition
      # compute them every {settings.check_termination} step
      mod(iter,settings.check_termination)  == 0 && ((r_prim,r_dual) = calculateResiduals(ws,settings))


      # check convergence with residuals every {settings.checkIteration} steps
      if mod(iter,settings.check_termination) == 0
        # update cost
        cost = ws.sm.cinv[]*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]

        if abs(cost) > 1e20
          status = :Unsolved
          break
        end

        # print iteration steps
        settings.verbose && printIteration(settings,iter,cost,r_prim,r_dual)

        if hasConverged(ws,settings,r_prim,r_dual)
          status = :Solved
          break
        end
      end

      # check infeasibility conditions every {settings.checkInfeasibility} steps
      if mod(iter,settings.check_infeasibility) == 0
        if isPrimalInfeasible(δy,ws,settings)
            status = :Primal_infeasible
            cost = Inf
            ws.x .= NaN
            ws.μ .= NaN
            ws.ν .= NaN
            break
        end

        if isDualInfeasible(δx,ws,settings)
            status = :Dual_infeasible
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

      if settings.time_limit !=0 &&  (time() - timeLimit_start) > settings.time_limit
        status = :Time_limit_reached
        break
      end

    end #END-ADMM-MAIN-LOOP

    settings.verbose_timing && (ws.times.iterTime = (time()-iter_start))
    settings.verbose_timing && (ws.times.postTime = time())

    # calculate primal and dual residuals
    if numIter == settings.max_iter
      r_prim,r_dual = calculateResiduals(ws,settings)
      status = :Max_iter_reached
    end

    # reverse scaling for scaled feasible cases
    if settings.scaling != 0 && (cost != Inf && cost != -Inf)
      reverseScaling!(ws)
      # FIXME: Another cost calculation is not necessary since cost value is not affected by scaling
      cost =  (1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1] #sm.cinv * not necessary anymore since reverseScaling
    end


    # print solution to screen

    ws.times.solverTime = time() - solverTime_start
    settings.verbose_timing && (ws.times.postTime = time()-ws.times.postTime)
    settings.verbose && printResult(status,numIter,cost,ws.times.solverTime)

    # create result object
    resinfo = QOCS.ResultInfo(r_prim,r_dual)
    y = -ws.μ

    return result = QOCS.Result(ws.x,y,ws.s,cost,numIter,status,resinfo,ws.times);


  end
