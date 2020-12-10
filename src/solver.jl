const LinsolveSubarray{T} = SubArray{T, 1, Vector{T},Tuple{UnitRange{Int}},true}

"""
	admm_z!
Evaluates the proximal operator of the Indicator function I_{R^n × K}.
"""
function admm_z!(x::Vector{T}, s::SplitVector{T},
	w::Vector{T},
	ρ::Vector{T},
	set::CompositeConvexSet{T}, m::Int64, n::Int64) where {T <: AbstractFloat}
	# 1) Projection step of w onto R^n x K
	@. x = w[1:n] #this is handled via a view, so no updating necessary

	@. s.data = w[n+1:end]
	p_time = @elapsed project!(s, set)

	# we recover μ from s and w
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
	mem = get_mem(ws.accelerator)
	# iter_history = COSMO.IterateHistory(m, n, mem)
	# COSMO.update_iterate_history!(iter_history, ws.vars.x, ws.vars.s, -ws.vars.μ, ws.vars.w, r_prim, r_dual, zeros(mem), NaN)
	

	COSMO.allocate_loop_variables!(ws, m, n)

	# warm starting the operator variable
	@. ws.vars.w[1:n] = ws.vars.x[1:n]
	@. ws.vars.w[n+1:n+m] = 1 / ws.ρvec * ws.vars.μ + ws.vars.s.data

	# change state of the workspace
	ws.states.IS_OPTIMIZED = true

	iter_start = time()

	# do one initialisation step to make ADMM iterates agree with standard ADMM
	COSMO.admm_x!(ws.vars.s, ws.ν, ws.s_tl, ws.ls, ws.sol, ws.vars.w, ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec, settings.sigma, m, n)
	COSMO.admm_w!(ws.vars.s, ws.x_tl, ws.s_tl, ws.vars.w, settings.alpha, m, n);

	for iter = 1:settings.max_iter
		num_iter += 1
	
		acceleration_pre!(ws.accelerator, ws, num_iter)
		
		if !was_succesful(ws.accelerator) && iter >= settings.check_infeasibility
			recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n) # μ_k kept in sync with s_k, w already updated to w_{k+1}
			@. ws.δy.data = ws.vars.μ
		end
			
		# ADMM steps
		@. ws.vars.w_prev = ws.vars.w
		ws.times.proj_time += admm_z!(ws.vars.x, ws.vars.s, ws.vars.w, ws.ρvec, ws.p.C, m, n) 
		admm_x!(ws.vars.s, ws.ν, ws.s_tl, ws.ls, ws.sol, ws.vars.w, ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec,settings.sigma, m, n)
		admm_w!(ws.vars.s, ws.x_tl, ws.s_tl, ws.vars.w, settings.alpha, m, n);	
		
		acceleration_post!(ws.accelerator, ws, num_iter)

		# 
		# check convergence with residuals every {settings.checkIteration} steps
		if mod(iter, settings.check_termination) == 0 || iter == 1
			recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)
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

		# computing δy requires one extra projection if acceleration is used, therefore update only
		# during times when acceleration wasn't used
		if !was_succesful(ws.accelerator) && iter > settings.check_infeasibility
			recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)
			@. ws.δy.data -= ws.vars.μ
		end
		# check infeasibility conditions every {settings.checkInfeasibility} steps
		if mod(iter, settings.check_infeasibility) == 0

			# compute δx for infeasibility detection
			@. ws.δx = ws.vars.w[1:n] - ws.vars.w_prev[1:n]

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

		# adapt ρvec at the appropriate intervals if 
		# - rho adaption is active {settings.adaptive_rho}
		# - rho has not been adapted {settings.adaptive_rho_max_adaptions} yet
		if settings.adaptive_rho && (settings.adaptive_rho_interval > 0) && (mod(iter, settings.adaptive_rho_interval) == 0) && (num_rho_adaptions(ws.rho_updates) < settings.adaptive_rho_max_adaptions)
			recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)
			was_adapted = adapt_rho_vec!(ws)
			# changing the rho changes the ADMM operator, so restart accelerator
			if was_adapted
				empty_history!(ws.accelerator)
				log_restart!(ws.accelerator, iter, :rho_adapted)
				# adapt w[n+1:end]
				@. ws.vars.w[n+1:end] = one(T) / ws.ρvec * ws.vars.μ + ws.vars.s.data
			end

		end

		if settings.time_limit !=0 &&  (time() - time_limit_start) > settings.time_limit
			status = :Time_limit_reached
			break
		end


	end #END-ADMM-MAIN-LOOP
	 
	recover_μ!(ws.vars.μ, ws.vars.w_prev, ws.vars.s, ws.ρvec, n)

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

	# if typeof(ws.accelerator) <: AndersonAccelerator{Float64}
	# # 	iter_history.aa_fail_data = ws.accelerator.fail_counter
	# end

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
	else
		@. ws.δx = zero(T)
		@. ws.δy.data = zero(T)
		@. ws.s_tl = zero(T)
		@. ws.ls = zero(T)
		@. ws.sol = zero(T)
	end

	ws.x_tl = view(ws.sol, 1:n) # i.e. xTilde
	ws.ν = view(ws.sol, (n + 1):(n + m))
end

"""
	acceleration_pre!(accelerator, ws, num_iter)

A function that the accelerator can use before the nominal ADMM operator step. 

In the case of an Anderson Accelerator this is used to calculate an accelerated candidate vector, that overwrites the current iterate.
"""
function acceleration_pre!(accelerator::AndersonAccelerator{T}, ws::Workspace{T}, num_iter::Int64) where {T <: AbstractFloat}
	COSMO.check_activation!(accelerator, num_iter)
	if is_active(accelerator)
		COSMO.update_history!(accelerator, ws.vars.w, ws.vars.w_prev, num_iter)
		# overwrite w here
		COSMO.accelerate!(ws.vars.w, ws.vars.w_prev, accelerator, num_iter)
	end 
end
acceleration_pre!(accelerator::AbstractAccelerator, args...; kwargs...) = nothing

"""
	acceleration_post!(accelerator, ws, num_iter)

A function that the accelerator can use after the nominal ADMM operator step. 

In the case of an Anderson Accelerator this is used to check the quality of the accelerated candidate vector and take measures if the vector is of bad quality.
"""
function acceleration_post!(accelerator::AndersonAccelerator{T}, ws::Workspace{T}, num_iter::Int64) where {T <: AbstractFloat}
	acc_post_time_start = time()
	if accelerator.success
		if is_safeguarding(accelerator) && is_active(accelerator)
			
			# norm(w_prev - w, 2) 
			# use accelerator.f here to get f = (x - g) as w_prev has been overwritten before this function call
			nrm_f = compute_non_accelerated_res_norm(accelerator.f)
			nrm_tol = nrm_f * accelerator.τ
			
			# compute residual norm of accelerated point
			nrm_f_acc = compute_accelerated_res_norm!(accelerator.f, ws.vars.w, ws.vars.w_prev)
			
			# Safeguarding check: f(w_acc) = w_prev_acc - w_acc <= τ f(w_prev) = τ * (w_prev - w) 
			# if safeguarding check fails, discard the accelerated candidate point and do a normal ADMM step instead
			if nrm_f_acc > nrm_tol
				accelerator.activate_logging && push!(accelerator.acceleration_status, (num_iter, :acc_guarded_declined))
				accelerator.num_safe_declined += 1
				# don't use accelerated candidate point. Reset w = g_last
				reset_accelerated_vector!(ws.vars.w, ws.vars.w_prev, accelerator.g_last)	
				m, n = ws.p.model_size 
				admm_z!(ws.vars.x, ws.vars.s, ws.vars.w, ws.ρvec, ws.p.C, m, n) 			
				admm_x!(ws.vars.s, ws.ν, ws.s_tl, ws.ls, ws.sol, ws.vars.w, ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec, ws.settings.sigma, m, n)
				admm_w!(ws.vars.s, ws.x_tl, ws.s_tl, ws.vars.w, ws.settings.alpha, m, n);

			else
				accelerator.num_accelerated_steps += 1
				accelerator.num_safe_accepted += 1
				accelerator.activate_logging && push!(accelerator.acceleration_status, (num_iter, :acc_guarded_accepted))

			end
			# For debugging: log the decision
			accelerator.activate_logging && push!(accelerator.safeguarding_status, (num_iter, nrm_f_acc, nrm_tol, nrm_f))
		else
			accelerator.activate_logging && push!(accelerator.acceleration_status, (num_iter, :acc_unguarded))
			accelerator.num_accelerated_steps += 1
		end
	end
	accelerator.acc_post_time += time() - acc_post_time_start
end

function compute_non_accelerated_res_norm(f::AbstractVector{T}) where {T <: AbstractFloat}
	return norm(f, 2)
end

function compute_accelerated_res_norm!(f::AbstractVector{T}, w::AbstractVector{T}, w_prev::AbstractVector{T}) where {T <: AbstractFloat}
	@. f = w_prev - w
	return norm(f, 2)	
end

function reset_accelerated_vector!(w::AbstractVector{T}, w_prev::AbstractVector{T}, g_last::AbstractVector{T}) where {T <: AbstractFloat}
		@. w_prev = g_last
		@. w = g_last
end


acceleration_post!(accelerator::AbstractAccelerator, args...; kwargs...) = nothing