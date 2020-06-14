const LinsolveSubarray{T} = SubArray{T, 1, Vector{T},Tuple{UnitRange{Int64}},true}


function admm_step!(x::Vector{T},
	s::SplitVector{T},
	μ::Vector{T},
	ν::LinsolveSubarray{T},
	x_tl::LinsolveSubarray{T},
	s_tl::Vector{T},
	ls::Vector{T},
	sol::Vector{T},
	kkt_solver::AbstractKKTSolver,
	q::Vector{T},
	b::Vector{T},
	ρ::Vector{T},
	α::T,
	σ::T,
	m::Int64,
	n::Int64,
	set::CompositeConvexSet{T}) where {T <: AbstractFloat}
	# linear solve
	# Create right hand side for linear system
	# deconstructed solution vector is ls = [x_tl(n+1); ν(n+1)]
	# x_tl and ν are automatically updated, since they are views on sol
	@. ls[1:n] = σ * x - q
	@. ls[(n + 1):end] = b - s.data + μ / ρ
	solve!(kkt_solver,sol,ls)

	# Over relaxation
	@. x = α * x_tl + (1.0 - α) * x
	@. s_tl = s.data - (ν + μ) / ρ
	@. s_tl = α * s_tl + (1.0 - α) * s.data
	@. s.data = s_tl + μ / ρ

	# Project onto cone
	p_time = @elapsed project!(s, set)

	# update dual variable μ
	@. μ = μ + ρ * (s_tl - s.data)
	return p_time
end

# SOLVER ROUTINE
# -------------------------------------


"""
optimize!(model)

Attempts to solve the optimization problem defined in `COSMO.Model` object with the user settings defined in `COSMO.Settings`. Returns a `COSMO.Result` object.
"""
function optimize!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}
	solver_time_start = time()
	settings = ws.settings
	# perform chordal decomposition
	if settings.decompose
		ws.times.graph_time = @elapsed COSMO.chordal_decomposition!(ws)
	end
	# create scaling variables
	# with scaling    -> uses mutable diagonal scaling matrices
	# without scaling -> uses identity matrices
	ws.sm = (settings.scaling > 0) ? ScaleMatrices{T}(ws.p.model_size[1], ws.p.model_size[2]) : ScaleMatrices{T}()

	# perform preprocessing steps (scaling, initial KKT factorization)

	# we measure times always in Float64
	ws.times.factor_update_time = 0.
	ws.times.proj_time  = 0. #reset projection time
	ws.times.setup_time = @elapsed setup!(ws);

	# instantiate variables
	status = :Unsolved
	cost = T(Inf)
	r_prim = T(Inf)
	r_dual = T(Inf)
	num_iter = 0
	# print information about settings to the screen
	settings.verbose && print_header(ws)
	time_limit_start = time()

	m, n = ws.p.model_size
	δx = zeros(T, n)
	δy = SplitVector{T}(zeros(T, m), ws.p.C)

	s_tl = zeros(T, m) # i.e. sTilde

	ls = zeros(T, n + m)
	sol = zeros(T, n + m)
	x_tl = view(sol, 1:n) # i.e. xTilde
	ν = view(sol, (n + 1):(n + m))

	iter_start = time()

	for iter = 1:settings.max_iter
		num_iter += 1

		# For infeasibility detection: Record the previous step just in time
		if mod(iter, settings.check_infeasibility) == settings.check_infeasibility - 1
			@. δx = ws.vars.x
			@. δy.data = ws.vars.μ
		end

		ws.times.proj_time += admm_step!(
			ws.vars.x, ws.vars.s, ws.vars.μ, ν,
			x_tl, s_tl, ls, sol,
			ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec,
			settings.alpha, settings.sigma,
			m, n, ws.p.C);

		# compute residuals (based on optimality conditions of the problem) to check for termination condition
		# compute them every {settings.check_termination} step


		# check convergence with residuals every {settings.checkIteration} steps
		if mod(iter, settings.check_termination) == 0 || iter == 1
			r_prim, r_dual = calculate_residuals!(ws)

			# update cost
			cost = calculate_cost!(ws.utility_vars.vec_n, ws.vars.x, ws.p.P, ws.p.q, ws.sm.cinv[])
			if abs(cost) > 1e20
				status = :Unsolved
				break
			end

			# print iteration steps
			settings.verbose && print_iteration(ws, iter, cost, r_prim, r_dual)

			if has_converged(ws, r_prim, r_dual)
				status = :Solved
				break
			end
		end

		# check infeasibility conditions every {settings.checkInfeasibility} steps
		if mod(iter, settings.check_infeasibility) == 0

			# compute deltas for infeasibility detection
			@. δx = ws.vars.x - δx
			@. δy.data -= ws.vars.μ

			if is_primal_infeasible!(δy, ws)
				status = :Primal_infeasible
				cost = Inf
				break
			end

			if is_dual_infeasible!(δx, ws)
				status = :Dual_infeasible
				cost = -Inf
				break
			end
		end

		# automatically choose adaptive rho interval if enabled in the settings
		if settings.adaptive_rho && settings.adaptive_rho_interval == 0
			# set the adaptive rho interval if a certain fraction of the setup time has passed
			if (time() - iter_start) > settings.adaptive_rho_fraction * ws.times.setup_time
				# round adaptive time interval to a multiple of the check_termination setting (and at least that)
				if settings.check_termination > 0
					settings.adaptive_rho_interval = round_multiple(iter, settings.check_termination)
					settings.adaptive_rho_interval = max(settings.adaptive_rho_interval, settings.check_termination)
				else
					settings.adaptive_rho_interval = round_multiple(iter, 25)
					settings.adaptive_rho_interval = max(settings.adaptive_rho_interval, 25)
				end
			end
		end

		# adapt rhoVec if enabled
		if settings.adaptive_rho && (settings.adaptive_rho_interval > 0) && (mod(iter, settings.adaptive_rho_interval) == 0)
			adapt_rho_vec!(ws)
		end

		if settings.time_limit !=0 &&  (time() - time_limit_start) > settings.time_limit
			status = :Time_limit_reached
			break
		end

	end #END-ADMM-MAIN-LOOP

	ws.times.iter_time = (time() - iter_start)
	settings.verbose_timing && (ws.times.post_time = time())

	# calculate primal and dual residuals
	if num_iter == settings.max_iter
		r_prim, r_dual = calculate_residuals!(ws)
		status = :Max_iter_reached
	end

	# reverse scaling for scaled feasible cases
	if settings.scaling != 0
		reverse_scaling!(ws)
		# FIXME: Another cost calculation is not necessary since cost value is not affected by scaling
		cost = calculate_cost!(ws.utility_vars.vec_n, ws.vars.x, ws.p.P, ws.p.q)
	end

	#reverse chordal decomposition
	if ws.ci.decompose
	 	reverse_decomposition!(ws, settings)
		y = -ws.vars.μ
 	else
		@. ws.utility_vars.vec_m = -ws.vars.μ
		y = ws.utility_vars.vec_m
	end

	ws.times.solver_time = time() - solver_time_start
	settings.verbose_timing && (ws.times.post_time = time() - ws.times.post_time)

	# print solution to screen
	settings.verbose && print_result(status, num_iter, cost, ws.times.solver_time)

	# create result object
	res_info = ResultInfo(r_prim, r_dual, ws.rho_updates)
	free_memory!(ws)

	return Result{T}(ws.vars.x, y, ws.vars.s.data, cost, num_iter, status, res_info, ws.times);

end


function free_memory!(ws)
	free_memory!(ws.kkt_solver)
end
