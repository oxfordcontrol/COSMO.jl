function factorKKT!(ws::COSMO.Workspace,settings::COSMO.Settings)

    p = ws.p

    # KKT matrix M
    M = [p.P+settings.sigma*I SparseMatrixCSC(p.A'); p.A -spdiagm(0=> (1.0./ws.ρVec))]
    # Do LDLT Factorization: A = LDL^T
    #try
    if settings.verbose_timing
        ws.times.factorTime += @elapsed ws.p.F = ldlt(M)
    else
        ws.p.F = ldlt(M)
    end
    #catch
    #  warn("Problems performing the LDLT facorization. Matrix has one or more zero pivots. Reusing previous step matrix.")
    #end
    return nothing
end
