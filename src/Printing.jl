module Printing
  using OSSDPTypes, Formatting
  export printHeader, printResult, printIteration

  function printHeader(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings,setupTime::Float64)
    n = ws.p.n
    m = ws.p.m
    K = ws.p.K
    numSOC = numSOCBlk = numSDP = numSDPBlk = 0
    if length(K.q) > 0
      numSOC = sum(K.q)
      numSOCBlk = length(K.q)
    else
      numSOC = 0
      numSOCBlk = 0
    end
    if length(K.s) > 0
      numSDP = sum(K.s)
      numSDPBlk = length(K.s)
    else
      numSDP = 0
      numSDPBlk = 0
    end
    settings.scaling > 0 ? scalingStatus = "on" : scalingStatus = "off"
    nnzInP = countnz(ws.p.P) - countnz(diag(ws.p.P)) + n
    nnzInM = 2*countnz(ws.p.A) + nnzInP + m

    println("-"^50 * "\n" * " "^8 * "ADMM-Conic Solver in pure Julia\n" * " "^18 * "Michael Garstka\n"  * " "^8 * "University of Oxford, February 2018\n" * "-"^50 * "\n")
    println("Problem:  x in R^{$(n)},\n          constraints: A in R^{$(n)x$(m)} ($(countnz(ws.p.A)) nnz), b in R^{$(m)},\n          matrix size to factor: $(n+m)x$(n+m) ($((n+m)^2) elem, $(nnzInM) nnz)\nCones:    free vars: $(ws.p.K.f)\n"*" "^10*"non-zero vars: $(ws.p.K.l)\n"*" "^10*"soc vars: $(numSOC), soc cones: $(numSOCBlk)\n"*" "^10*"sdp vars: $(numSDP), sdp cones: $(numSDPBlk)")
    println("Settings: ϵ_abs = $(settings.eps_abs), ϵ_rel = $(settings.eps_rel),\n" * " "^10 * "ϵ_prim_inf = $(settings.eps_prim_inf), ϵ_dual_inf = $(settings.eps_dual_inf),\n" * " "^10 * "ρ = $(settings.rho), σ = $(settings.sigma), α = $(settings.alpha),\n" * " "^10 * "max_iter = $(settings.max_iter), scaling = $(settings.scaling) ($(scalingStatus))")
    println("Setup Time: $(round.(setupTime*1000,2))ms\n")
  end

  function printIteration(settings::OSSDPTypes.OSSDPSettings,iter::Int64,cost::Float64,r_prim::Float64,r_dual::Float64)
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


  function printResult(status::Symbol,iter::Int64,cost::Float64,rt::Float64)
    println("\n" * "-"^50 * "\nRESULT: Status: $(status)\nTotal Iterations: $(iter)\nOptimal objective: $(round.(cost,4))\nRuntime: $(round.(rt,3))s ($(round.(rt*1000,2))ms)\n" * "-"^50 )
  end
end