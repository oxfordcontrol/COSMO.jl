using Printf

function print_header(ws::COSMO.Workspace)

	m, n = ws.p.model_size
	settings = ws.settings
	settings.scaling > 0 ? scaling_status = "on" : scaling_status = "off"
	nnz_in_P = count(!iszero,ws.p.P) - count(!iszero,diag(ws.p.P)) + n
	nnz_in_M = 2*count(!iszero,ws.p.A) + nnz_in_P + m
	println("-"^66 * "\n" * " "^10 * "COSMO v0.6.0 - A Quadratic Objective Conic Solver\n" * " "^25 * "Michael Garstka\n"  * " "^16 * "University of Oxford, 2017 - 2019\n" * "-"^66 * "\n")

	println("Problem:  x ∈ R^{$(n)},\n          constraints: A ∈ R^{$(m)x$(n)} ($(count(!iszero, ws.p.A)) nnz),\n          matrix size to factor: $(n + m)x$(n + m) ($(nnz_in_M) nnz)")
	for (iii, set) in enumerate(sort(ws.p.C.sets, by = x -> -x.dim))
		set_name = split(string(typeof(set).name), ".")[end]
		if isa(set, PsdCone) || isa(set, PsdConeTriangle) || isa(set, PsdConeTriangleLOBPCG)
			set_dim = string(set.sqrt_dim, "x", set.sqrt_dim)
		else
			set_dim = string(set.dim)
		end
		iii == 1 ? println("Sets:"*" "^5*"$(set_name) of dim: $(set_dim)") : println(" "^10*"$(set_name) of dim: $(set_dim)")
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
	if settings.psd_projector == PsdConeTriangleLOBPCG || (isa(settings.psd_projector, OptionsFactory) && settings.psd_projector.ObjectType == PsdConeTriangleLOBPCG)
		println("Iter:\tObjective:\tPrimal Res:\tDual Res:\tRho:\t\tTime(s):\tSubspaces:\tLOBPCG avg it:\tExact avg it:")
	else
		println("Iter:\tObjective:\tPrimal Res:\tDual Res:\tRho:\t\tTime(s):")
	end
	flush(stdout)
	nothing
end

function print_iteration(ws::COSMO.Workspace, iter::Int64, cost::Float64, r_prim::Float64, r_dual::Float64, running_time::Float64)
	settings = ws.settings
	if mod(iter, 1) == 0 || iter == 1 || iter == 2 || iter == settings.max_iter
		if mod(iter, settings.check_termination) == 0
			lobpcg_dims_str, lobpcg_iterations, exact_projections = extract_lobpcg_statistics(ws)
			if isnan(lobpcg_iterations)
				@printf("%d\t%.4e\t%.4e\t%.4e\t%.4e\t%.2e\n", iter, cost, r_prim, r_dual, ws.ρ, running_time)
			else
				@printf("%d\t%.4e\t%.4e\t%.4e\t%.4e\t%.2e\t%s\t%.2f\t\t%.2f\n", iter, cost, r_prim, r_dual, ws.ρ, running_time,
					lobpcg_dims_str, lobpcg_iterations, exact_projections)
			end
		else
			@printf("%d\t%.4e\t ---\t\t\t---\t%.2e\n", iter, cost, running_time)
		end
	end
	flush(stdout)
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

function extract_lobpcg_statistics(ws; max_printed_dimensions=3)
	lobpcg_iterations, exact_projections = 0, 0
	lobpcg_dims = Int[]
	sets_counter = 0
	for set in ws.p.C.sets
		if isa(set, COSMO.PsdConeTriangleLOBPCG)
			if length(lobpcg_dims) < max_printed_dimensions
				lobpcg_dim = set.subspace_dim_history[end]
				if set.is_subspace_positive
					append!(lobpcg_dims, lobpcg_dim)
				else
					append!(lobpcg_dims, -lobpcg_dim)
				end
			end

			lobpcg_iterations += set.lobpcg_iterations
			exact_projections += set.exact_projections
			set.lobpcg_iterations = 0
			set.exact_projections = 0
			sets_counter += 1
		end
	end
	lobpcg_dims_str = string(lobpcg_dims)
	if length(lobpcg_dims) < sets_counter
		lobpcg_dims_str = string(lobpcg_dims_str[1:end-1], ", ...]")
	end
	lobpcg_dims_str = rpad(lobpcg_dims_str, 12, " ")
	return lobpcg_dims_str, lobpcg_iterations/sets_counter, exact_projections/sets_counter
end


function print_result(status::Symbol, iter::Int64, cost::Float64, rt::Float64)
	println("\n" * "-"^66 * "\n>>> Results\nStatus: $(status)\nIterations: $(iter)\nOptimal objective: $(round.(cost; digits = 4))\nRuntime: $(round.(rt; digits = 3))s ($(round.(rt * 1000; digits = 2))ms)\n")
	nothing
end

function print_lin_sys(s::Type{<:AbstractKKTSolver})
	lin_sys_dict = Dict("QdldlKKTSolver" => "QDLDL", "CholmodKKTSolver" => "CHOLMOD", "PardisoDirectKKTSolver" => "Pardiso (direct)", "PardisoIndirectKKTSolver" => "Pardiso (indirect)",  "MKLPardisoKKTSolver" => "MKL Pardiso")
	return get(lin_sys_dict, string(s), string(s))
end