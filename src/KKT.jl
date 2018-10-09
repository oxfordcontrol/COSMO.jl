module KKT
using ..QOCS, LinearAlgebra, SparseArrays
export factorKKT!

  function factorKKT!(ws::QOCS.Workspace,settings::QOCS.Settings)
     p = ws.p
     if nnz(p.P) > 0 && p.P != p.P'
      p.P = p.P ./ 2+(p.P ./ 2)'
    end
    # KKT matrix M
    M = [p.P+settings.sigma*sparse(1.0I,p.n,p.n) p.A';p.A -sparse(Diagonal(1 ./ws.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    #try
    if settings.verbose_timing
      ws.times.factorTime += @elapsed ws.p.F = ldlt(M)
      @show ws.times.factorTime
    else
      ws.p.F = ldlt(M)
    end
    #catch
    #  warn("Problems performing the LDLT facorization. Matrix has one or more zero pivots. Reusing previous step matrix.")
    #end
    return nothing
  end

end #MODULE