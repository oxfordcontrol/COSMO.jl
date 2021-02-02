const LinsolveSubarray{T} = SubArray{T, 1, Vector{T},Tuple{UnitRange{Int}},true}

"""
	admm_z!
Evaluates the proximal operator of the Indicator function I_{R^n × K}.
"""
function admm_z!(s::SplitVector{T},
	w::Vector{T},
	set::CompositeConvexSet{T}, n::Int64) where {T <: AbstractFloat}
	# 1) Projection step of w onto R^n x K
	# @. x = w[1:n] #this is handled via a view: x = view(w_prev, 1:n), prev to keep it in sync with s

	# 2) s = Π(w)
	@. s.data = w[n+1:end]
	p_time = @elapsed project!(s, set)

	# 3) y = ρ * (w - Π(w)) (Moreau decomposition)
	# we recover μ from s and w just-in-time
	# @. μ = ρ * (w[n+1:end] - s.data)
	return p_time
end

"The dual variable μ can be recovered from w, s via Moreau decomposition: μ = ρ (w - Π(w))."
function recover_μ!(μ::Vector{T}, w::Vector{T}, s::SplitVector{T}, ρ::Vector{T}, n::Int64) where {T <: AbstractFloat}
	@. μ = ρ * (w[n+1:end] - s.data)
end

"""
	admm_x!
Evaluates the proximal operator of 1/2 x'Px + q'x s.t. Ax + s == b.
"""
function admm_x!(s::SplitVector{T},
	ν::LinsolveSubarray{T},
	s_tl::Vector{T},
	ls::Vector{T},
	sol::Vector{T},
	w::Vector{T},
	kkt_solver::AbstractKKTSolver,
	q::Vector{T},
	b::Vector{T},
	ρ::Vector{T},	
	σ::T,
	m::Int,
	n::Int) where {T <: AbstractFloat}

	# 2) linear solve
	# Create right hand side for linear system
	# deconstructed solution vector is ls = [x_tl(n+1); ν(n+1)]
	# x_tl and ν are automatically updated, since they are views on sol
	@. ls[1:n] = σ * w[1:n] - q	# T(2) * σ * x - σ * w[1:n] - q, as x == w[1:n] at that point
	@. ls[(n + 1):end] = b - T(2) * s.data + w[(n + 1):end]
	solve!(kkt_solver, sol, ls)

	# x_tl and ν are automatically updated as they are views into sol
	@. s_tl = T(2) * s.data - w[n+1:end] - ν  / ρ
end

"""
	admm_w!
ADMM-operator variable `w` update with over-relaxation parameter α.
"""
function admm_w!(s::SplitVector{T}, x_tl::LinsolveSubarray{T}, s_tl::Vector{T}, w::Vector{T}, α::T, m::Int64, n::Int64) where {T <: AbstractFloat}
	@. w[1:n] = w[1:n] + α * (x_tl - w[1:n]) 
	@. w[n+1:end] = w[n+1:end] + α * (s_tl - s.data)
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
	if settings.verbose_timing
		ws.times.update_time = 0.
		ws.times.accelerate_time = 0.
	end
	
	# instantiate variables
	status = :Undetermined
	cost = T(Inf)
	r_prim = T(Inf)
	r_dual = T(Inf)
	iter = 0
	ws.safeguarding_iter = 0
	# print information about settings to the screen
	settings.verbose && print_header(ws)
	time_limit_start = time()

	m, n = ws.p.model_size


	COSMO.allocate_loop_variables!(ws, m, n)

	# warm starting the operator variable
	@. ws.vars.w[1:n] = ws.vars.x[1:n]
	@. ws.vars.w[n+1:n+m] = one(T) / ws.ρvec * ws.vars.μ + ws.vars.s.data

	# change state of the workspace
	ws.states.IS_OPTIMIZED = true

	iter_start = time()

	# do one initialisation step to make ADMM iterates agree with standard ADMM
	COSMO.admm_x!(ws.vars.s, ws.ν, ws.s_tl, ws.ls, ws.sol, ws.vars.w, ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec, settings.sigma, m, n)
	COSMO.admm_w!(ws.vars.s, ws.x_tl, ws.s_tl, ws.vars.w, settings.alpha, m, n);

	while iter + ws.safeguarding_iter < settings.max_iter
		iter += 1

		acceleration_pre!(ws.accelerator, ws, iter)
		
		if update_suggested(ws.infeasibility_check_due, ws.accelerator)
			recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n) # μ_k kept in sync with s_k, w already updated to w_{k+1}
			@. ws.δy.data = ws.vars.μ
		end

		# ADMM steps
		@. ws.vars.w_prev = ws.vars.w
		ws.times.proj_time += admm_z!(ws.vars.s, ws.vars.w, ws.p.C, n) 
		apply_rho_adaptation_rules!(ws.ρvec, ws.rho_updates, settings, iter, iter_start, ws.times, ws, n)
		admm_x!(ws.vars.s, ws.ν, ws.s_tl, ws.ls, ws.sol, ws.vars.w, ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec,settings.sigma, m, n)
		admm_w!(ws.vars.s, ws.x_tl, ws.s_tl, ws.vars.w, settings.alpha, m, n);	

		acceleration_post!(ws.accelerator, ws, iter)

		# convergence / infeasibility / timelimit checks
		cost, status, r_prim, r_dual = check_termination!(ws, settings, iter, cost, status, r_prim, r_dual, time_limit_start, n)
		if status != :Undetermined
			break
		end

	end #END-ADMM-MAIN-LOOP
	 
	recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)

	ws.times.iter_time = (time() - iter_start)
	settings.verbose_timing && (ws.times.post_time = time())

	# calculate primal and dual residuals
	if iter + ws.safeguarding_iter == settings.max_iter
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
	total_iter = ws.safeguarding_iter + iter
	settings.verbose && print_result(status, total_iter, ws.safeguarding_iter, cost, ws.times.solver_time, ws.settings.safeguard)

	# create result object
	res_info = ResultInfo(r_prim, r_dual, ws.rho_updates)
	free_memory!(ws)
	return Result{T}(ws.vars.x, y, ws.vars.s.data, cost, total_iter, ws.safeguarding_iter, status, res_info, ws.times);

end

"Free all memory before exiting optimize!"
function free_memory!(ws)
	free_memory!(ws.kkt_solver)
end

"Allocate ADMM- helper variables once the problem dimension is determined."
function allocate_loop_variables!(ws::COSMO.Model{T}, m::Int, n::Int) where {T <: AbstractFloat}

	if length(ws.δx) != n || length(ws.ls) != m + n
		ws.δx = zeros(T, n)
		ws.δy = SplitVector{T}(zeros(T, m), ws.p.C)
		ws.s_tl = zeros(T, m)
		ws.ls = zeros(T, n + m)
		ws.sol = zeros(T, n + m)
	else
		@. ws.δx = zero(T)
		@. ws.δy.data = zero(T)
		@. ws.s_tl = zero(T)
		@. ws.ls = zero(T)
		@. ws.sol = zero(T)
	end

	ws.x_tl = view(ws.sol, 1:n) # i.e. xTilde
	ws.ν = view(ws.sol, (n + 1):(n + m))
	ws.rho_update_due = false
end

"""
	apply_rho_adaptation_rules!()

Adapts the step-size parameter ρ when multiple conditions are satisfied.

- If automatic rho interval is active, the `adaptive_rho_interval` is chosen as a fraction of the setup time.

- Adapt rho if the iteration count reaches a multiple of `adaptive_rho_interval` and the maximum number of allowable updates has not been reached (`adaptive_rho_max_adaptions`).
- If an `accelerator` is used, update rho at the next possible iterations, i.e. at the next non-accelerated iteration.
"""
function apply_rho_adaptation_rules!(ρvec::Vector{T}, rho_updates::Vector{T}, settings::Settings{T}, iter::Int64, iter_start::Float64, times::ResultTimes{Float64}, ws::Workspace{T}, n::Int64) where {T <: AbstractFloat}
	# automatically choose adaptive rho interval if enabled in the settings
	if settings.adaptive_rho && settings.adaptive_rho_interval == 0
		# set the adaptive rho interval if a certain fraction of the setup time has passed
		if (time() - iter_start) > settings.adaptive_rho_fraction * times.setup_time
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

	# adapt ρvec at the appropriate intervals if 
	# - rho adaption is active {settings.adaptive_rho}
	# - rho has not been adapted {settings.adaptive_rho_max_adaptions} times yet
	if settings.adaptive_rho && (settings.adaptive_rho_interval > 0) && (mod(iter, settings.adaptive_rho_interval) == 0) && (num_rho_adaptions(rho_updates) < settings.adaptive_rho_max_adaptions)
		ws.rho_update_due = true
	end
	

	# adapt rho at the next possible iteration
	# when an accelerator is used this will take place at the next non-accelerated iteration
	if update_suggested(ws.rho_update_due, ws.accelerator)
		ws.rho_update_due = false
		recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ρvec, n)
		was_adapted = adapt_rho_vec!(ws)
		# changing the rho changes the ADMM operator, so restart accelerator
		if was_adapted
			CA.restart!(ws.accelerator)
			CA.log!(ws.accelerator, iter, :rho_adapted)
			
			# adapt w[n+1:end]
			@. ws.vars.w[n+1:end] = one(T) / ρvec * ws.vars.μ + ws.vars.s.data
		end

	end
end

update_suggested(update_due::Bool, aa::AbstractAccelerator) = update_due

function update_suggested(update_due::Bool, aa::AndersonAccelerator)
 		if update_due && !was_successful(aa)
		return true
	else
		return false
	end
end

"""
	check_termination!()

Checks the algorithms termination conditions at intervals specified in the `settings` and updates the algorithm `status`:
- residuals are "small enough" --> :Solved
- cost value out-ouf-range --> :Unsolved
- infeasibility conditions satisfied --> :Primal_infeasible / :Dual_infeasible
- time limit constraint reached --> :Time_limit_reached
"""
function check_termination!(ws::Workspace{T}, settings::Settings{T}, iter::Int64, cost::T, status::Symbol, r_prim::T, r_dual::T, time_limit_start::Float64, n::Int64) where {T <: AbstractFloat}

	# check convergence with residuals every {settings.checkIteration} steps
	if mod(iter, settings.check_termination) == 0 || iter == 1
		recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)
		r_prim, r_dual = calculate_residuals!(ws)
		# update cost
		cost = calculate_cost!(ws.utility_vars.vec_n, ws.vars.x, ws.p.P, ws.p.q, ws.sm.cinv[])
		if abs(cost) > 1e20
			status = :Unsolved
			return cost, status, r_prim, r_dual
		end
		# print iteration steps
		settings.verbose && print_iteration(ws, iter, cost, r_prim, r_dual)
		if has_converged(ws, r_prim, r_dual)
			status = :Solved
			return cost, status, r_prim, r_dual	
		end
	end

	# computing δy requires one extra projection if acceleration is used, therefore update only
	# during times when acceleration wasn't used, i.e. set the ws.infeasibility_check_due flag
	# and δy, δx will be computed and used for the check at the next non-accelerated iteration 
	if mod(iter, settings.check_infeasibility) == 0
		ws.infeasibility_check_due = true
	else
		if update_suggested(ws.infeasibility_check_due, ws.accelerator)
			ws.infeasibility_check_due = false
			recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)
			@. ws.δy.data -= ws.vars.μ

			# compute δx for infeasibility detection
			@. ws.δx = ws.vars.w[1:n] - ws.vars.w_prev[1:n]

			if is_primal_infeasible!(ws.δy, ws)
				status = :Primal_infeasible
				cost = Inf
				return cost, status, r_prim, r_dual
			end

			if is_dual_infeasible!(ws.δx, ws)
				status = :Dual_infeasible
				cost = -Inf
				return cost, status, r_prim, r_dual
			end
		end
	end

	if settings.time_limit !=0 &&  (time() - time_limit_start) > settings.time_limit
		status = :Time_limit_reached
	end
	return cost, status, r_prim, r_dual
end