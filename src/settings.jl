export with_options

struct OptionsFactory{T}
	# See https://github.com/JuliaOpt/JuMP.jl/blob/94e2cbb3b1a079b327d2fb9b80bfd3b3afa3415d/src/JuMP.jl#L84-L121
	ObjectType::Type
	args::Tuple
	kwargs::Base.Iterators.Pairs
end

function with_options(object_type, args...; kwargs...)
	OptionsFactory{object_type}(object_type, args, kwargs)
end

function (options_factory::OptionsFactory{T})(args...; kwargs...) where {T}
    return options_factory.ObjectType(args..., options_factory.args...;
									kwargs..., options_factory.kwargs...)
end

"""
	COSMO.Settings(; kws)

Creates a COSMO settings object that is used to pass user settings to the solver.

Argument | Description | Values (default)
--- | :--- | :---
rho | ADMM rho step | 0.1
sigma | ADMM sigma step | 1e-6
alpha | Relaxation parameter | 1.6
eps_abs | Absolute residual tolerance | 1e-4
eps_rel | Relative residual tolerance | 1e-4
eps\\_prim\\_inf | Primal infeasibility tolerance | 1e-4
eps\\_dual\\_inf | Dual infeasibility tolerance | 1e-4
max_iter | Maximum number of iterations | 2500
verbose | Verbose printing | false
verbose_timing | Verbose timing | false
kkt_solver | Linear System solver | `CholmodKKTSolver`
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
decompose | Activate to decompose chordal psd constraints | true
complete_dual | Activate to complete the dual variable after decomposition | false
merge_strategy | Choose a strategy for clique merging | `CliqueGraphMerge`
compact_transformation | Choose how a decomposed problem is transformed | true
time_limit | Set solver time limit in s | 0
"""
mutable struct Settings{T <: AbstractFloat}
	rho::T
	sigma::T
	alpha::T
	eps_abs::T
	eps_rel::T
	eps_prim_inf::T
	eps_dual_inf::T
	max_iter::Int64
	verbose::Bool
	kkt_solver::Union{Type{<: AbstractKKTSolver}, OptionsFactory{<: AbstractKKTSolver}}
	check_termination::Int64
	check_infeasibility::Int64
	scaling::Int64
	MIN_SCALING::T
	MAX_SCALING::T
	adaptive_rho::Bool
	adaptive_rho_interval::Int64
	adaptive_rho_tolerance::T
	adaptive_rho_fraction::T
	verbose_timing::Bool
	RHO_MIN::T
	RHO_MAX::T
	RHO_TOL::T
	RHO_EQ_OVER_RHO_INEQ::T
	COSMO_INFTY::T
  	decompose::Bool
  	complete_dual::Bool
	time_limit::T
	obj_true::T
	obj_true_tol::T
	merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}}
	compact_transformation::Bool
	#constructor
	function Settings{T}(;
		rho=T(0.1),
		sigma=T(1e-6),
		alpha=T(1.6),
		eps_abs=T(1e-4),
		eps_rel=T(1e-4),
		eps_prim_inf=T(1e-6),
		eps_dual_inf=T(1e-4),
		max_iter=2500,
		verbose=false,
		kkt_solver=QdldlKKTSolver,
		check_termination=40,
		check_infeasibility=40,
		scaling=10,
		MIN_SCALING = T(1e-4),
		MAX_SCALING = T(1e4),
		adaptive_rho = true,
		adaptive_rho_interval = 40,
		adaptive_rho_tolerance = 5,
		adaptive_rho_fraction = T(0.4),
		verbose_timing = false,
		RHO_MIN = T(1e-6),
		RHO_MAX = T(1e6),
		RHO_TOL = T(1e-4),
		RHO_EQ_OVER_RHO_INEQ = T(1e3),
		COSMO_INFTY = T(1e20),
		decompose = true,
    	complete_dual = false,
		time_limit = 0.0,
		obj_true = NaN,
		obj_true_tol = T(1e-3),
		merge_strategy = CliqueGraphMerge,
		compact_transformation = true
		) where {T <: AbstractFloat}
		if !isa(kkt_solver, OptionsFactory)
			kkt_solver = with_options(kkt_solver)
		end

		if !isa(merge_strategy, OptionsFactory)
			merge_strategy = with_options(merge_strategy)
		end
		new(rho, sigma, alpha, eps_abs, eps_rel, eps_prim_inf, eps_dual_inf, max_iter, verbose, kkt_solver, check_termination, check_infeasibility, scaling, MIN_SCALING, MAX_SCALING, adaptive_rho, adaptive_rho_interval, adaptive_rho_tolerance, adaptive_rho_fraction, verbose_timing, RHO_MIN, RHO_MAX, RHO_TOL, RHO_EQ_OVER_RHO_INEQ, COSMO_INFTY, decompose, complete_dual, time_limit, obj_true, obj_true_tol, merge_strategy, compact_transformation)
	end
end

Settings(args...) = Settings{DefaultFloat}(args...)
