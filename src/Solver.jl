include("./Types.jl")
include("./Scaling.jl")
include("./Projections.jl")
include("./Residuals.jl")
include("./Parameters.jl")
include("./Infeasibility.jl")
include("./Printing.jl")
include("./KKT.jl")
include("./Setup.jl")

module OSSDP

using Formatting, Projections, Scaling, OSSDPTypes, Parameters, Infeasibility, Residuals, Printing, Setup
export solve, OSSDPSettings, Cone #from the Types module


# SOLVER ROUTINE
# -------------------------------------
  function solve(P,q,A,b,K::OSSDPTypes.Cone,settings::OSSDPTypes.OSSDPSettings)

    # create workspace variables
    ws = WorkSpace(Problem(P,q,A,b,K),ScaleMatrices())
    P = q = A = b = nothing

    # perform preprocessing steps (scaling, initial KKT factorization)
    setup!(ws,settings)

    # instantiate variables
    iter = 0
    status = :unsolved
    cost = Inf
    δν = []

    # print information about settings to the screen
    settings.verbose && printHeader(ws,settings)

    tic()

    # MAIN ADMM LOOP
    for iter = 1:1:settings.max_iter

      # assign previous variables, here x(n+1) becomes the new x(n), s(n+1) -> x(n)

      # construct right hand side [-q+σ*s(n)-λ(n); b-(1/ρ)*ν(n)]
      RHS = [-ws.p.q+settings.sigma*ws.s-ws.λ; ws.p.b-diagm((1./ws.p.ρVec))*ws.ν]

      #solve linear system M*k = b with help of factorization matrix
      # FIXME: must be a better way
      k = sparse(ws.p.F\full(RHS))

      # The relaxation definitely has to be double checked
      #deconstruct solution vector k = [x(n+1);ν(n+1)]
      ws.x = k[1:ws.p.n]
      νPrev = ws.ν
      ws.ν = k[ws.p.n+1:end]

      # Projection steps and relaxation xRelax = αx(n+1)+(1-α)s(n)
      xRelax = settings.alpha*ws.x+(1-settings.alpha)*ws.s

      # s(n+1) = Proj( xRelax + (1/σ)*λ(n))
      ws.s = Projections.projectCompositeCone!((xRelax + (1/settings.sigma)*ws.λ),ws.p.K)

      # update dual variables λ(n+1) = λ(n) + σ*(xRelax - s(n+1))
      ws.λ = ws.λ + settings.sigma*(xRelax - ws.s)

      # update cost
      # FIXME: Remove calculation of cost at each step
      cost = ws.sm.cinv*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]

      # compute residuals (based on optimality conditions of the problem) to check for termination condition
      # compute them every {settings.checkTermination} step
      mod(iter,settings.checkTermination)  == 0 && ((r_prim,r_dual) = calculateResiduals(ws,settings))

      # compute deltas
      # δx = xNew - xPrev
      # δy = yNew - yPrev
      push!(δν,norm(ws.ν-νPrev,Inf))
      # print iteration steps
      if settings.verbose
        if iter == 1
          println("Iter:\tObjective:\tPrimal Res\tDual Res:")
        end
        if mod(iter,1) == 0 || iter == 1 || iter == 2 || iter == settings.max_iter
          if mod(iter,settings.checkTermination) == 0
            printfmt("{1:d}\t{2:.4e}\t{3:.4e}\t{4:.4e}\n", iter,cost,r_prim,r_dual)
          else
            printfmt("{1:d}\t{2:.4e}\t ---\t\t\t---\n", iter,cost)
          end
       end
      end

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
          if settings.verbose
            printfmt("{1:d}\t{2:.4e}\t{3:.4e}\t{4:.4e}\n", iter,cost,r_prim,r_dual)
          end
          status = :solved
          break
        end
      end

      # adapt rhoVec if enabled
      if settings.adaptive_rho && (mod(iter,settings.adaptive_rho_interval) == 0) && (settings.adaptive_rho_interval > 0)
        adaptRhoVec!(ws,settings)
      end

    end #END-ADMM-MAIN-LOOP

    rt = toq()

    # calculate primal and dual residuals
    if iter == settings.max_iter
      r_prim,r_dual = calculateResiduals(ws,settings)
      status = :UserLimit
    end

    xS = copy(ws.x)
    sS = copy(ws.s)
    νS = copy(ws.ν)
    λS = copy(ws.λ)


    if settings.scaling != 0
      reverseScaling!(ws)
      # FIXME: Another cost calculation is not necessary since cost value is not affected by scaling
      cost =  (1/2 * ws.x'*p.P*ws.x + p.q'*ws.x)[1] #sm.cinv * not necessary anymore since reverseScaling
    end

    # print solution to screen
    settings.verbose && printResult(status,iter,cost,rt)


    # create result object
    result = OSSDPResult(ws.x,ws.s,ws.λ,ws.ν,cost,iter,status,rt,r_prim,r_dual);

    return result,ws,xS,sS,νS,λS,δν;

  end

end




