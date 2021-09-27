# functions that interact with the accelerator
export IterActivation, AccuracyActivation, ImmediateActivation
const CA = COSMOAccelerators

"Activate accelerator after `start_iter` iterations of the main algorithm."
struct IterActivation <: AbstractActivationReason
  start_iter::Int64
  function IterActivation(start_iter::Int64)
    start_iter < 2 && ArgumentError("start_iter has to be at least 2.")
    new(start_iter)
  end
end

"Activate accelerator after accuracy of main algorithm <= `start_accuracy`."
struct AccuracyActivation <: AbstractActivationReason
  start_accuracy::Float64
  function AccuracyActivation(start_accuracy::Float64)
    start_accuracy >= 0. && ArgumentError("start_accuracy has to be a non-negative number.")
    new(start_accuracy)
  end 
end

# dispatch on the activation reason
function check_activation!(ws::Workspace, activation_reason::ImmediateActivation, num_iter::Int64)
    if !ws.accelerator_active && num_iter >= 2
        ws.accelerator_active = true
    end
end

function check_activation!(ws::Workspace, activation_reason::IterActivation, num_iter::Int64) 
    if !ws.accelerator_active && num_iter >= activation_reason.start_iter
        ws.accelerator_active = true
    end
end
check_activation!(ws::Workspace, activation_reason::AccuracyActivation, num_iter::Int64) = nothing


function check_activation!(ws::Workspace, activation_reason::AccuracyActivation, r::ResultInfo) where {T <: AbstractFloat}
    if !ws.accelerator_active
        tol = activation_reason.start_accuracy
        
        if isapprox_primal_feasible(r, tol, tol) && isapprox_dual_feasible(r, tol, tol)
            ws.accelerator_active = true
        end 
    end
end

check_activation!(ws::Workspace, activation_reason::Union{IterActivation, ImmediateActivation}, r::ResultInfo) where {T <: AbstractFloat}= nothing


"""
	acceleration_pre!(accelerator, ws, num_iter)

A function that the accelerator can use before the nominal ADMM operator step. 

In the case of an Anderson Accelerator this is used to calculate an accelerated candidate vector, that overwrites the current iterate.
"""
function acceleration_pre!(accelerator::CA.AbstractAccelerator, ws::Workspace{T}, num_iter::Int64) where {T <: AbstractFloat}
	COSMO.check_activation!(ws, ws.activation_reason, num_iter)
	if ws.accelerator_active
		if ws.settings.verbose_timing
			time_start = time()  
			CA.update!(accelerator, ws.vars.w, ws.vars.w_prev, num_iter)
			ws.times.update_time += time() - time_start
			time_start = time()  		
			CA.accelerate!(ws.vars.w, ws.vars.w_prev, accelerator, num_iter)
			ws.times.accelerate_time += time() - time_start
		else
			CA.update!(accelerator, ws.vars.w, ws.vars.w_prev, num_iter)
			# overwrite w here
			CA.accelerate!(ws.vars.w, ws.vars.w_prev, accelerator, num_iter)
		end
	end 
end
acceleration_pre!(accelerator::CA.EmptyAccelerator, ws::Workspace{T}, num_iter::Int64) where {T <: AbstractFloat} = nothing


"""
	acceleration_post!(accelerator, ws, num_iter, safeguarding_iter)

A function that the accelerator can use after the nominal ADMM operator step. 

In the case of an Anderson Accelerator this is used to check the quality of the accelerated candidate vector and take measures if the vector is of bad quality.
"""
function acceleration_post!(accelerator::CA.AndersonAccelerator, ws::Workspace{T}, num_iter::Int64) where {T <: AbstractFloat}
	if ws.accelerator_active && CA.was_successful(accelerator)
		if ws.settings.safeguard
			# norm(w_prev - w, 2) 
			# use accelerator.f here to get f = (x - g) as w_prev has been overwritten before this function call
			nrm_f = norm(accelerator.f, 2)
			nrm_tol = nrm_f * ws.settings.safeguard_tol
			
			# compute residual norm of accelerated point
			nrm_f_acc = compute_accelerated_res_norm!(accelerator.f, ws.vars.w, ws.vars.w_prev)
			
			# Safeguarding check: f(w_acc) = w_prev_acc - w_acc <= τ f(w_prev) = τ * (w_prev - w) 
			# if safeguarding check fails, discard the accelerated candidate point and do a normal ADMM step instead
            if nrm_f_acc > nrm_tol
                CA.log!(accelerator, num_iter, :acc_guarded_declined)
				# don't use accelerated candidate point. Reset w = g_last
				reset_accelerated_vector!(ws.vars.w, ws.vars.w_prev, accelerator)	
				m, n = ws.p.model_size 
				ws.times.proj_time += admm_z!(ws.vars.s, ws.vars.w, ws.p.C, n) 			
				admm_x!(ws.vars.s, ws.ν, ws.s_tl, ws.ls, ws.sol, ws.vars.w, ws.kkt_solver, ws.p.q, ws.p.b, ws.ρvec, ws.settings.sigma, m, n)
				admm_w!(ws.vars.s, ws.x_tl, ws.s_tl, ws.vars.w, ws.settings.alpha, m, n);
				ws.safeguarding_iter += 1
            else
                CA.log!(accelerator, num_iter, :acc_guarded_accepted)
			end
        else
            CA.log!(accelerator, num_iter, :acc_unguarded) 
		end
	end
end

acceleration_post!(accelerator::CA.AbstractAccelerator, ws::Workspace{T}, num_iter::Int64) where {T <: AbstractFloat} = nothing


"Computes the nonexpansive operator residual norm, i.e. f(w_prev) = || w_prev - w ||_2."
function compute_accelerated_res_norm!(f::AbstractVector{T}, w::AbstractVector{T}, w_prev::AbstractVector{T}) where {T <: AbstractFloat}
	@. f = w_prev - w
	return norm(f, 2)	
end

"Reset ADMM variables to last non-accelerated point."
reset_accelerated_vector!(w::AbstractVector{T}, w_prev::AbstractVector{T}, aa::CA.AndersonAccelerator{T}) where {T <: AbstractFloat} = reset_accelerated_vector!(w, w_prev, aa.g_last)
function reset_accelerated_vector!(w::AbstractVector{T}, w_prev::AbstractVector{T}, g_last::AbstractVector{T}) where {T <: AbstractFloat}
		@. w_prev = g_last
		@. w = g_last
end
