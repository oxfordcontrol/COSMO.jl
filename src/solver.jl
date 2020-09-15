const LinsolveSubarray{T} = SubArray{T, 1, Vector{T},Tuple{UnitRange{Int}},true}


function admm_step!(x::Vector{T},
	s::SplitVector{T},
	μ::Vector{T},
	ν::LinsolveSubarray{T},
	x_tl::LinsolveSubarray{T},
	s_tl::Vector{T},
	ls::Vector{T},
	sol::Vector{T},
	w::Vector{T},
	kkt_solver::AbstractKKTSolver,
	q::Vector{T},
	b::Vector{T},
	ρ::Vector{T},
	α::T,
	σ::T,
	m::Int,
	n::Int,
	set::CompositeConvexSet{T}) where {T <: AbstractFloat}

	# 1) Projection step of w
	@. x = w[1:n]
	@. s = w[n+1:end]
	# Project onto cone
	p_time = @elapsed project!(s, set)
	# we recover μ from s and w
	@. μ = ρ * (w[n+1:end] - s)

	# 2) linear solve
	# Create right hand side for linear system
	# deconstructed solution vector is ls = [x_tl(n+1); ν(n+1)]
	# x_tl and ν are automatically updated, since they are views on sol
	@. ls[1:n] = T(2) * σ * x - σ * w[1:n] - q
	@. ls[(n + 1):end] = b - T(2) * s.data + w[(n + 1):end]
	solve!(kkt_solver, sol, ls)

	# x_tl and ν are automatically updated as they are views into sol
	@. s_tl = T(2) * s.data - w[n+1:end] - ν  / ρ

	# 3) dual variable update with over-relaxation
	@. w[1:n] = w[1:n] + α * (x_tl - x)
	@. w[n+1:end] = w[n+1:end] + α * (s_tl - s)
	return p_time
end

# SOLVER ROUTINE
# -------------------------------------


"""
optimize!(model)

Attempts to solve the optimization problem defined in `COSMO.Model` object with the user settings defined in `COSMO.Settings`. Returns a `COSMO.Result` object.
"""
function optimize!(ws::COSMO.Workspace{T}) where {T <: AbstractFloat}

	!ws.states.IS_ASSEMBLED && throw(ErrorException("The model has to be assembled! / set! before optimize!() can be called."))

	# start timer
	solver_time_start = time()

	settings = ws.settings

	# perform chordal decomposition
	if settings.decompose
		if !ws.states.IS_CHORDAL_DECOMPOSED
			ws.times.graph_time = @elapsed COSMO.chordal_decomposition!(ws)
		elseif ws.ci.decompose
			pre_allocate_variables!(ws)
		end
	end

	# create scaling variables
	# with scaling    -> uses mutable diagonal scaling matrices
	# without scaling -> uses identity matrices
	if !ws.states.IS_SCALED
		ws.sm = (settings.scaling > 0) ? COSMO.ScaleMatrices{T}(ws.p.model_size[1], ws.p.model_size[2]) : COSMO.ScaleMatrices{T}()
	end

	# we measure times always in Float64
	ws.times.factor_update_time = 0.
	ws.times.proj_time  = 0. #reset projection time
	ws.times.setup_time = @elapsed COSMO.setup!(ws);

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
	iter_history = IterateHistory(m, n, settings.acc_mem)
	update_iterate_history!(iter_history, ws.vars.x, ws.vars.s, -ws.vars.μ, ws.vars.w, r_prim, r_dual, zeros(ws.settings.acc_mem), NaN)

	COSMO.allocate_loop_variables!(ws, m, n)

	@. w[1:n] = ws.vars.x[1:n]
	@. w[n+1:end] = ws.vars.s.data

	x_tl = view(sol, 1:n) # i.e. xTilde
	ν = view(sol, (n + 1):(n + m))
	iter_start = time()

	for iter = 1:settings.max_iter
		num_iter += 1

		if num_iter > 1
			update_history!(ws.accelerator, ws.vars.w, ws.vars.w_prev)
			accelerate!(ws.vars.w, ws.vars.w_prev, ws.accelerator, num_iter)
		end

		# For infeasibility detection: Record the previous step just in time
		if mod(iter, settings.check_infeasibility) == settings.check_infeasibility - 1
			@. ws.δx = ws.vars.x
			@. ws.δy.data = ws.vars.μ
		end

		ws.times.proj_time += COSMO.admm_step!(
			ws.vars.x, ws.vars.s, ws.vars.μ, ν,
			x_tl, ws.s_tl, ws.ls, ws.sol, ws.w,
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
			@. ws.δx = ws.vars.x - ws.δx
			@. ws.δy.data -= ws.vars.μ

			if is_primal_infeasible!(ws.δy, ws)
				status = :Primal_infeasible
				cost = Inf
				break
			end

			if is_dual_infeasible!(ws.δx, ws)
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

	if typeof(ws.accelerator) <: AndersonAccelerator{Float64}
		iter_history.aa_fail_data = ws.accelerator.fail_counter
	end

	return Result{T}(ws.vars.x, y, ws.vars.s.data, cost, num_iter, status, res_info, ws.times);

end


function free_memory!(ws)
	free_memory!(ws.kkt_solver)
end

function allocate_loop_variables!(ws::COSMO.Model{T}, m::Int, n::Int ) where {T <: AbstractFloat}

	if length(ws.δx) != n || length(ws.ls) != m + n
		ws.δx = zeros(T, n)
		ws.δy = SplitVector{T}(zeros(T, m), ws.p.C)
		ws.s_tl = zeros(T, m)
		ws.ls = zeros(T, n + m)
		ws.sol = zeros(T, n + m)
		ws.w = zeros(T, n + m) # the α-averaged variable of the algorithm
	else
		@. ws.δx = zero(T)
		@. ws.δy.data = zero(T)
		@. ws.s_tl = zero(T)
		@. ws.ls = zero(T)
		@. ws.sol = zero(T)
		@. ws.w = zero(T)
	end

end
