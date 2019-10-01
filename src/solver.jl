const LinsolveSubarray = SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true}

function admm_step!(x::Vector{Float64},
	s::SplitVector{Float64},
	μ::Vector{Float64},
	ν::LinsolveSubarray,
	x_tl::LinsolveSubarray,
	s_tl::Vector{Float64},
	ls::Vector{Float64},
	sol::Vector{Float64},
	kkt_solver::AbstractKKTSolver,
	q::Vector{Float64},
	b::Vector{Float64},
	ρ::Vector{Float64},
	α::Float64,
	σ::Float64,
	m::Int64,
	n::Int64,
	set::CompositeConvexSet{Float64})
	# linear solve
	# Create right hand side for linear system
	# deconstructed solution vector is ls = [x_tl(n+1); ν(n+1)]
	# x_tl and ν are automatically updated, since they are views on sol
	@. ls[1:n] = σ * x - q
	@. ls[(n + 1):end] = b - s + μ / ρ
	solve!(kkt_solver,sol,ls)

	# Over relaxation
	@. x = α * x_tl + (1.0 - α) * x
	@. s_tl = s - (ν + μ) / ρ
	@. s_tl = α * s_tl + (1.0 - α) * s
	@. s = s_tl + μ / ρ

	# Project onto cone
	p_time = @elapsed project!(s, set)

	# update dual variable μ
	@. μ = μ + ρ .* (s_tl - s)
	return p_time
end

# SOLVER ROUTINE
# -------------------------------------


"""
optimize!(model)

Attempts to solve the optimization problem defined in `COSMO.Model` object with the user settings defined in `COSMO.Settings`. Returns a `COSMO.Result` object.
"""
function optimize!(ws::COSMO.Workspace)
	solver_time_start = time()
	settings = ws.settings
  # perform chordal decomposition
  if settings.decompose
  	ws.times.graph_time = @elapsed COSMO.chordal_decomposition!(ws)
  end
	# create scaling variables
	# with scaling    -> uses mutable diagonal scaling matrices
	# without scaling -> uses identity matrices
	ws.sm = (settings.scaling > 0) ? ScaleMatrices(ws.p.model_size[1], ws.p.model_size[2]) : ScaleMatrices()

	# perform preprocessing steps (scaling, initial KKT factorization)
	ws.times.factor_time = 0
	ws.times.proj_time  = 0. #reset projection time
	ws.times.setup_time = @elapsed setup!(ws);

	# instantiate variables
	num_iter = 0
	status = :Unsolved
	cost = Inf
	r_prim = Inf
	r_dual = Inf

	# print information about settings to the screen
	settings.verbose && print_header(ws)
	time_limit_start = time()

	m, n = ws.p.model_size
	δx = zeros(n)
	δy = SplitVector(zeros(m), ws.p.C)

	s_tl = zeros(m) # i.e. sTilde

	ls = zeros(n + m)
	sol = zeros(n + m)
	x_tl = view(sol, 1:n) # i.e. xTilde
	ν = view(sol, (n + 1):(n + m))

	settings.verbose_timing && (iter_start = time())

	for iter = 1:settings.max_iter

		num_iter+= 1

		@. δx = ws.vars.x
		@. δy = ws.vars.μ

		ws.times.proj_time += admm_step!(
			ws.vars.x, ws.vars.s, ws.vars.μ, ν,
			x_tl, s_tl, ls, sol,
			ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec,
			settings.alpha, settings.sigma,
			m, n, ws.p.C);

		# compute deltas for infeasibility detection
		@. δx = ws.vars.x - δx

		#@. δy = -ws.vars.μ + δy
		@. δy -= ws.vars.μ


		# compute residuals (based on optimality conditions of the problem) to check for termination condition
		# compute them every {settings.check_termination} step
		mod(iter, settings.check_termination)  == 0 && ((r_prim, r_dual) = calculate_residuals!(ws))

		# check convergence with residuals every {settings.checkIteration} steps
		if mod(iter, settings.check_termination) == 0
			# update cost
			cost = ws.sm.cinv[] * (1/2 * ws.vars.x' * ws.p.P * ws.vars.x + ws.p.q' * ws.vars.x)[1]

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

		# adapt rhoVec if enabled
		if settings.adaptive_rho && (mod(iter, settings.adaptive_rho_interval) == 0) && (settings.adaptive_rho_interval > 0)
			adapt_rho_vec!(ws)
		end

		if settings.time_limit !=0 &&  (time() - time_limit_start) > settings.time_limit
			status = :Time_limit_reached
			break
		end

	end #END-ADMM-MAIN-LOOP

	settings.verbose_timing && (ws.times.iter_time = (time() - iter_start))
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
		cost =  (1/2 * ws.vars.x' * ws.p.P * ws.vars.x + ws.p.q' * ws.vars.x)[1] #sm.cinv * not necessary anymore since reverseScaling
	end

	#reverse chordal decomposition
	if settings.decompose
	 reverse_decomposition!(ws, settings)
	end

	ws.times.solver_time = time() - solver_time_start
	settings.verbose_timing && (ws.times.post_time = time() - ws.times.post_time)
	# print solution to screen
	settings.verbose && print_result(status, num_iter, cost, ws.times.solver_time)

	# create result object
	res_info = ResultInfo(r_prim, r_dual)

	y = -ws.vars.μ
	free_memory!(ws)

	return Result{Float64}(ws.vars.x, y, ws.vars.s.data, cost, num_iter, status, res_info, ws.times), ws;

end


function free_memory!(ws)
	free_memory!(ws.kkt_solver)
end
