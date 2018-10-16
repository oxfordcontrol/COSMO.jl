function factorKKT!(ws::COSMO.Workspace,settings::COSMO.Settings)

    p = ws.p
    @show "Malakia"

    if length(p.M) == 0
        # KKT matrix
        p.M = [p.P+settings.sigma*I SparseMatrixCSC(p.A'); p.A -I]
    end
    n = p.model_size[2]
    m = p.model_size[1]
    @inbounds @simd for i = (n + 1):(n + m)
        p.M[i, i] = -1.0/ws.ρVec[i - n]
    end

    # Do LDLT Factorization: A = LDL^T
    #try
    if settings.verbose_timing
        ws.times.factorTime += @elapsed ws.p.F = ldlt(p.M)
    else
        ws.p.F = ldlt(p.M)
    end
    #catch
    #  warn("Problems performing the LDLT facorization. Matrix has one or more zero pivots. Reusing previous step matrix.")
    #end
    return nothing
end
