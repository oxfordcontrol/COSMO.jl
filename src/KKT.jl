module KKT
using OSSDPTypes, Helper
export factorKKT!

  function factorKKT!(p::OSSDPTypes.Problem,settings::OSSDPTypes.OSSDPSettings)
     if nnz(p.P) > 0 && p.P != p.P'
      i,j,difference = findNonSymmetricComponent(p.P)
      warn("Scaled P is not symmetric. [$(i),$(j)] differs by $(difference). Trying to correct.")
      p.P = p.P./2+(p.P./2)'
    end
    # KKT matrix M
    M = [p.P+settings.sigma*speye(p.n) p.A';p.A -spdiagm((1./p.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    p.F = ldltfact(M)
    # p.F = M
    return nothing
  end

end #MODULE