module KKT
using ..QOCS, LinearAlgebra, SparseArrays
export factorKKT!

  function factorKKT!(ws::QOCS.Workspace,settings::QOCS.Settings)
     p = ws.p
    if length(p.M) <= 1
      At = SparseMatrixCSC(p.A')
      p.M = [p.P+settings.sigma*I At; p.A -I]
    end
    @inbounds @simd for i = (p.n + 1):(p.n + p.m)
      p.M[i, i] = -1.0/ws.ρVec[i - ws.p.n]
    end
    # Do LDLT Factorization: A = LDL^T
    #try
    if settings.verbose_timing
      ws.times.factorTime += @elapsed ws.p.F = ldlt(p.M)
      @show ws.times.factorTime
    else
      ws.p.F = ldlt(p.M)
    end
    #catch
    #  warn("Problems performing the LDLT facorization. Matrix has one or more zero pivots. Reusing previous step matrix.")
    #end
    return nothing
  end

end #MODULE