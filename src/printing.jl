print_version() = println(COSMO.version())

function print_header(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}

	m, n = ws.p.model_size
	settings = ws.settings
	settings.scaling > 0 ? scaling_status = "on" : scaling_status = "off"

	println("-"^66 * "\n" * " "^10 * "COSMO v0.8.9 - A Quadratic Objective Conic Solver\n" * " "^25 * "Michael Garstka\n"  * " "^16 * "University of Oxford, 2017 - 2022\n" * "-"^66 * "\n")

	println("Problem:  x ∈ R^{$(n)},\n          constraints: A ∈ R^{$(m)x$(n)} ($(count(!iszero, ws.p.A)) nnz),\n          matrix size to factor: $(n + m)x$(n + m),\n          Floating-point precision: $(T)")
	for (iii, set) in enumerate(sort(ws.p.C.sets, by = x -> -x.dim))
		if iii == 1
			print("Sets:"*" "^5)
			print_set(set)
		else
			print(" "^10)
			print_set(set)
		end
		if iii == 5
			println(" "^10*"... and $(length(ws.p.C.sets)-5) more")
			break
		end
	end
	if ws.ci.decompose && ws.ci.num_decomposable > 0
		println("Decomp:   Num of original PSD cones: $(ws.ci.num_psd_cones)\n" * " "^10*"Num decomposable PSD cones: $(ws.ci.num_decomposable)\n" * " "^10*"Num PSD cones after decomposition: $(ws.ci.num_decom_psd_cones)\n" * " "^10*"Merge Strategy: $(stringify(ws.settings.merge_strategy))")
	end
	println("Settings: ϵ_abs = $(@sprintf("%.1e", settings.eps_abs)), ϵ_rel = $(@sprintf("%.1e", settings.eps_rel)),\n" * " "^10 * "ϵ_prim_inf = $(@sprintf("%.1e", settings.eps_prim_inf)), ϵ_dual_inf = $(@sprintf("%.1e", settings.eps_dual_inf)),\n" * " "^10 * "ρ = $(@sprintf("%.3g", settings.rho)), σ = $(@sprintf("%.3g", settings.sigma)), α = $(@sprintf("%.3g", settings.alpha)),\n" * " "^10 * "max_iter = $(settings.max_iter),\n" * " "^10 * "scaling iter = $(settings.scaling) ($(scaling_status)),\n" * " "^10 * "check termination every $(settings.check_termination) iter,\n" * " "^10 * "check infeasibility every $(settings.check_infeasibility) iter,\n" *	 " "^10 * "KKT system solver: $(print_lin_sys(settings.kkt_solver))")

	print_accelerator(ws.accelerator, ws.settings.safeguard, ws.settings.safeguard_tol, tab = 10)
	

	println("Setup Time: $(round.(ws.times.setup_time*1000; digits=2))ms\n")
	println("Iter:\tObjective:\tPrimal Res:\tDual Res:\tRho:")
	nothing
end

function print_iteration(ws::COSMO.Workspace, iter::Int, cost::T, r_prim::T, r_dual::T) where {T <: AbstractFloat}
	settings = ws.settings
	if mod(iter, 1) == 0 || iter == 1 || iter == 2 || iter == settings.max_iter
		if mod(iter, settings.check_termination) == 0 || iter == 1
			@printf("%d\t% .4e\t%.4e\t%.4e\t%.4e\n", iter, cost, r_prim, r_dual, ws.ρ)
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

function print_result(status::Symbol, iter::Int, safeguarding_iter::Int, cost::T, rt::Float64, safeguard::Bool) where {T <: AbstractFloat}
	print("\n" * "-"^66 * "\n>>> Results\nStatus: ")
	if status == :Solved
		result_color = :green
	else
		result_color = :red
	end
	printstyled("$(status)\n", color = result_color)
	print("Iterations: $(iter)")
	if safeguard && safeguarding_iter > 0
		println(" (incl. $(safeguarding_iter) safeguarding iter)")
	else
		print("\n")
	end

	println("Optimal objective: $(@sprintf("%.4g", cost))\nRuntime: $(round.(rt; digits = 3))s ($(round.(rt * 1000; digits = 2))ms)\n")
	nothing
end
print_lin_sys(s::OptionsFactory) = print_lin_sys(s.ObjectType)
function print_lin_sys(s::Type{<:AbstractKKTSolver})
	lin_sys_dict = Dict("QdldlKKTSolver" => "QDLDL", "CholmodKKTSolver" => "CHOLMOD", "PardisoDirectKKTSolver" => "Pardiso (direct)", "PardisoIndirectKKTSolver" => "Pardiso (indirect)",  "MKLPardisoKKTSolver" => "MKL Pardiso")
	return get(lin_sys_dict, string(s), string(s))
end

"Print information about the accelerator at the start."
function print_accelerator(s::CA.AndersonAccelerator{T, BT, ME, RE}, safeguard::Bool, safeguard_tol; tab::Int64) where {T, BT, ME, RE}
		me = ME == CA.RestartedMemory ? "RestartedMemory" : "RollingMemory"
		if BT == CA.Type1
			bt = "Type1"
		elseif BT == CA.Type2{NormalEquations}
			bt = "Type2{NormalEquations}"
		elseif BT == CA.Type2{QRDecomp}
			bt = "Type2{QRDecomp}"	
		end

		println("Acc:" * " "^(tab - 4)  * "Anderson $(bt),\n" * " "^tab * "Memory size = $(s.mem), $(me),	\n" * " "^tab * "Safeguarded: $(safeguard), tol: $(safeguard_tol)")
end


print_accelerator(s::CA.AbstractAccelerator, safeguarded::Bool, safeguarding_tol; tab::Int64) = println("Acc:" * " "^(tab - 4) * "Unknown Accelerator")
print_accelerator(s::CA.EmptyAccelerator,  safeguarded::Bool, safeguarding_tol; tab::Int64) = println("Acc:" * " "^(tab - 4) * "no acceleration")

function print_set(set::COSMO.AbstractConvexSet)

	set_name = split(string(typeof(set).name), ".")[end][1:end-1]
	println("$(set_name) of dim: $(set.dim)")
end

function print_set(set::Union{PsdCone{T}, DensePsdCone{T}, PsdConeTriangle{T}, DensePsdConeTriangle{T}}) where {T <: AbstractFloat}
	set_name = split(string(typeof(set).name), ".")[end][1:end-1]
	N = set.sqrt_dim
	println("$(set_name) of dim: $(set.dim) ($(N)x$(N))")
end
