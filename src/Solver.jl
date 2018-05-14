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

# SOLVER ROUTINE
# -------------------------------------
  function solve(P,q,A,b,K::OSSDPTypes.Cone,settings::OSSDPTypes.OSSDPSettings)

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

    tic()
    startTime = time()
    xtilde = zeros(ws.p.n)
    stilde = zeros(ws.p.m)

    # MAIN ADMM LOOP
    for iter = 1:1:settings.max_iter

      # assign previous variables, here x(n+1) becomes the new x(n), s(n+1) -> s(n)

      # construct right hand side [-q+σ*x(n); b-s(n)+(1/ρ)μ(n)]
      RHS = [-ws.p.q+settings.sigma*ws.x; ws.p.b-ws.s+diagm(1./ws.p.ρVec)*ws.μ]

      #solve linear system M*k = b with help of factorization matrix
      # FIXME: must be a better way
      k = sparse(ws.p.F\full(RHS))

      #deconstruct solution vector k = [x(n+1);ν(n+1)]
      xtilde = k[1:ws.p.n]
      ws.ν = k[ws.p.n+1:end]
      stilde = ws.s - diagm(1./ws.p.ρVec)*(ws.ν+ws.μ)

      # Projection steps
      ws.x = settings.alpha*xtilde + (1-settings.alpha)*ws.x
      stildeRelax = settings.alpha*stilde + (1-settings.alpha)*ws.s
      ws.s = Projections.projectCompositeCone!(stildeRelax + diagm(1./ws.p.ρVec)*ws.μ,ws.p.K)

      # update dual variable μ
      ws.μ = ws.μ + diagm(ws.p.ρVec)*(stildeRelax - ws.s)

      # update cost
      # FIXME: Remove calculation of cost at each step
      cost = ws.sm.cinv*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]

      # compute residuals (based on optimality conditions of the problem) to check for termination condition
      # compute them every {settings.checkTermination} step
      mod(iter,settings.checkTermination)  == 0 && ((r_prim,r_dual) = calculateResiduals(ws,settings))

      # compute deltas
      # δx = xNew - xPrev
      # δy = yNew - yPrev

      # print iteration steps
      settings.verbose && printIteration(settings,iter,cost,r_prim,r_dual)


      # if isPrimalInfeasible(δy,A,l,u,settings.ϵ_prim_inf)
      #     status = :primal_infeasible
      #     cost = Inf
      #     xNew = NaN*ones(n,1)
      #     yNew = NaN*ones(m,1)
      #     break
      # end

      # if isDualInfeasible(δx,P,A,q,l,u,settings.ϵ_dual_inf)
      #     status = :dual_infeasible
      #     cost = -Inf
      #     xNew = NaN*ones(n,1)
      #     yNew = NaN*ones(m,1)
      #     break
      # end

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




