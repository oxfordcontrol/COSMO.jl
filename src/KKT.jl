module KKT
using OSSDPTypes
export factorKKT!

  function factorKKT!(p::OSSDPTypes.Problem,settings::OSSDPTypes.OSSDPSettings)
     if p.P != p.P'
      warn("Scaled P is not symmetric. Trying to correct.")
      p.P = p.P./2+(p.P./2)'
    end
    # KKT matrix M
    M = [p.P+settings.sigma*speye(p.n) p.A';p.A -spdiagm((1./p.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    p.F = ldltfact(M)
    return nothing
  end

end #MODULE