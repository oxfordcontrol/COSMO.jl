module KKT
using QOCS
export factorKKT!

  function factorKKT!(p::QOCS.Problem,settings::QOCS.Settings)
     if nnz(p.P) > 0 && p.P != p.P'
      p.P = p.P./2+(p.P./2)'
    end
    # KKT matrix M
    M = [p.P+settings.sigma*speye(p.n) p.A';p.A -spdiagm((1./p.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    #try
      p.F = ldltfact(M)
    #catch
    #  warn("Problems performing the LDLT facorization. Matrix has one or more zero pivots. Reusing previous step matrix.")
    #end
    return nothing
  end

end #MODULE