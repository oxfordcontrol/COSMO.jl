module KKT
using OSSDPTypes
export factorKKT!

  function factorKKT!(p::OSSDPTypes.Problem,settings::OSSDPTypes.OSSDPSettings)
    # KKT matrix M
    M = [p.P+settings.sigma*speye(p.n) p.A';p.A -spdiagm((1./p.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    p.F = ldltfact(M)
    return nothing
  end

end #MODULE