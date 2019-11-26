using Printf

function print_header(ws::COSMO.Workspace)

	m, n = ws.p.model_size
	settings = ws.settings
	settings.scaling > 0 ? scaling_status = "on" : scaling_status = "off"
	nnz_in_P = count(!iszero,ws.p.P) - count(!iszero,diag(ws.p.P)) + n
	nnz_in_M = 2*count(!iszero,ws.p.A) + nnz_in_P + m
	println("-"^66 * "\n" * " "^10 * "COSMO v0.5.0 - A Quadratic Objective Conic Solver\n" * " "^25 * "Michael Garstka\n"  * " "^16 * "University of Oxford, 2017 - 2019\n" * "-"^66 * "\n")

	println("Problem:  x ∈ R^{$(n)},\n          constraints: A ∈ R^{$(m)x$(n)} ($(count(!iszero, ws.p.A)) nnz),\n          matrix size to factor: $(n + m)x$(n + m) ($(nnz_in_M) nnz)")
	for (iii, set) in enumerate(sort(ws.p.C.sets, by = x -> -x.dim))
		set_name = split(string(typeof(set).name), ".")[end]
		iii == 1 ? println("Sets:"*" "^5*"$(set_name) of dim: $(set.dim)") : println(" "^10*"$(set_name) of dim: $(set.dim)")
		if iii > 5
			println(" "^10*"... and $(length(ws.p.C.sets)-5) more")
			break
		end
	end
	if ws.settings.decompose && ws.ci.num_decomposable > 0
		println("Decomp:   Num of original PSD cones: $(ws.ci.num_psd_cones)\n" * " "^10*"Num decomposable PSD cones: $(ws.ci.num_decomposable)\n" * " "^10*"Num PSD cones after decomposition: $(ws.ci.num_decom_psd_cones)\n" * " "^10*"Merge Strategy: $(stringify(ws.settings.merge_strategy))")
	end
	println("Settings: ϵ_abs = $(@sprintf("%.1e", settings.eps_abs)), ϵ_rel = $(@sprintf("%.1e", settings.eps_rel)),\n" * " "^10 * "ϵ_prim_inf = $(@sprintf("%.1e", settings.eps_prim_inf)), ϵ_dual_inf = $(@sprintf("%.1e", settings.eps_dual_inf)),\n" * " "^10 * "ρ = $(settings.rho), σ = $(settings.sigma), α = $(settings.alpha),\n" * " "^10 * "max_iter = $(settings.max_iter),\n" * " "^10 * "scaling iter = $(settings.scaling) ($(scaling_status)),\n" * " "^10 * "check termination every $(settings.check_termination) iter,\n" * " "^10 * "check infeasibility every $(settings.check_infeasibility) iter,\n" *	 " "^10 * "KKT system solver: $(print_lin_sys(settings.kkt_solver.ObjectType))")

	println("Setup Time: $(round.(ws.times.setup_time*1000; digits=2))ms\n")
	println("Iter:\tObjective:\tPrimal Res:\tDual Res:\tRho:")
	nothing
end

function print_iteration(ws::COSMO.Workspace, iter::Int64, cost::Float64, r_prim::Float64, r_dual::Float64)
	settings = ws.settings
	if mod(iter, 1) == 0 || iter == 1 || iter == 2 || iter == settings.max_iter
		if mod(iter, settings.check_termination) == 0
			@printf("%d\t%.4e\t%.4e\t%.4e\t%.4e\n", iter, cost, r_prim, r_dual, ws.ρ)
		else
			@printf("%d\t%.4e\t ---\t\t\t---\n", iter, cost)
		end
	end
	nothing
end

function stringify(merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}})
	if typeof(merge_strategy) <: OptionsFactory
		return stringify(merge_strategy.ObjectType)
	else
		s = string(merge_strategy)
		# remove the COSMO prefix
		return s[7:end]
	end
end

function print_result(status::Symbol, iter::Int64, cost::Float64, rt::Float64)
	println("\n" * "-"^66 * "\n>>> Results\nStatus: $(status)\nIterations: $(iter)\nOptimal objective: $(round.(cost; digits = 4))\nRuntime: $(round.(rt; digits = 3))s ($(round.(rt * 1000; digits = 2))ms)\n")
	nothing
end

function print_lin_sys(s::Type{<:AbstractKKTSolver})
	lin_sys_dict = Dict("QdldlKKTSolver" => "QDLDL", "CholmodKKTSolver" => "CHOLMOD", "PardisoDirectKKTSolver" => "Pardiso (direct)", "PardisoIndirectKKTSolver" => "Pardiso (indirect)",  "MKLPardisoKKTSolver" => "MKL Pardiso")
	return get(lin_sys_dict, string(s), string(s))
end