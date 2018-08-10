module KKT
using ..QOCS, LinearAlgebra, SparseArrays
export factorKKT!

  function factorKKT!(p::QOCS.Problem,settings::QOCS.Settings)
     if nnz(p.P) > 0 && p.P != p.P'
      p.P = p.P ./ 2+(p.P ./ 2)'
    end
    # KKT matrix M
    M = [p.P+settings.sigma*sparse(1.0I,p.n,p.n) p.A';p.A -sparse(Diagonal(1 ./p.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    #try
      p.F = ldlt(M)
    #catch
    #  warn("Problems performing the LDLT facorization. Matrix has one or more zero pivots. Reusing previous step matrix.")
    #end
    return nothing
  end

end #MODULE