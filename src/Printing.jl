module Printing
  using OSSDPTypes
  export printHeader, printResult

  function printHeader(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings)
    n = ws.p.n
    m = ws.p.m
    println("-"^50 * "\n" * " "^8 * "ADMM-SDP Solver in pure Julia\n" * " "^18 * "Michael Garstka\n"  * " "^8 * "University of Oxford, February 2018\n" * "-"^50 * "\n")
    println("Problem: variable X in ???, vec(X) in R^{$(n)},\n         constraints: A in R^{$(n)x$(m)}, b in R^{$(m)},\n         matrix size to factor: $(n+m)x$(n+m) ($((n+m)^2) elem)")
    println("Settings: ϵ_abs = $(settings.eps_abs), ϵ_rel = $(settings.eps_rel),\n" * " "^10 * "ϵ_prim_inf = $(settings.eps_prim_inf), ϵ_dual_inf = $(settings.eps_dual_inf),\n" * " "^10 * "ρ = $(settings.rho), σ = $(settings.sigma), α = $(settings.alpha),\n" * " "^10 * "max_iter = $(settings.max_iter), scaling = $(settings.scaling)\n\n")
  end


  function printResult(status,iter,cost,rt)
    println("\n\n" * "-"^50 * "\nRESULT: Status: $(status)\nTotal Iterations: $(iter)\nOptimal objective: $(round.(cost,4))\nRuntime: $(round.(rt,3))s ($(round.(rt*1000,2))ms)\n" * "-"^50 )
  end
end