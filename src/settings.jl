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
	Settings(;arg=val)

Creates a COSMO settings object that is used to pass user settings to the solver.

Argument | Description | Values (default)
--- | --- | ---
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
kkt_solver | Linear System solver | QDLDLKKTSolver
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
decompose | Activate to decompose chordal psd constraints | false
complete_dual | Activate to complete the dual variable after decomposition | false
merge_strategy | Choose a strategy for clique merging | PairwiseMerge
time_limit | set solver time limit in s | 0
"""
mutable struct Settings
	rho::Float64
	sigma::Float64
	alpha::Float64
	eps_abs::Float64
	eps_rel::Float64
	eps_prim_inf::Float64
	eps_dual_inf::Float64
	max_iter::Int64
	verbose::Bool
	kkt_solver::Union{Type{<: AbstractKKTSolver}, OptionsFactory{<: AbstractKKTSolver}}
	check_termination::Int64
	check_infeasibility::Int64
	scaling::Int64
	MIN_SCALING::Float64
	MAX_SCALING::Float64
	adaptive_rho::Bool
	adaptive_rho_interval::Int64
	adaptive_rho_tolerance::Float64
	verbose_timing::Bool
	RHO_MIN::Float64
	RHO_MAX::Float64
	RHO_TOL::Float64
  decompose::Bool
  complete_dual::Bool
	time_limit::Float64
	obj_true::Float64
	obj_true_tol::Float64
	merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}}
	colo_transformation::Bool
	#constructor
	function Settings(;
		rho=0.1,
		sigma=1e-6,
		alpha=1.6,
		eps_abs=1e-4,
		eps_rel=1e-4,
		eps_prim_inf=1e-6,
		eps_dual_inf=1e-4,
		max_iter=2500,
		verbose=false,
		kkt_solver=QdldlKKTSolver,
		check_termination=40,
		check_infeasibility=40,
		scaling=10,
		MIN_SCALING = 1e-4,
		MAX_SCALING = 1e4,
		adaptive_rho = true,
		adaptive_rho_interval = 40,
		adaptive_rho_tolerance = 5,
		verbose_timing = false,
		RHO_MIN = 1e-6,
		RHO_MAX = 1e6,
		RHO_TOL = 1e-4,
		decompose = false,
    complete_dual = false,
		time_limit = 0.0,
		obj_true = NaN,
		obj_true_tol = 1e-3,
		merge_strategy = PairwiseMerge,
		colo_transformation = true
		)
	if !isa(kkt_solver, OptionsFactory)
		kkt_solver = with_options(kkt_solver)
	end

	if !isa(merge_strategy, OptionsFactory)
		merge_strategy = with_options(merge_strategy)
	end
	new(rho, sigma, alpha, eps_abs, eps_rel, eps_prim_inf, eps_dual_inf, max_iter, verbose, kkt_solver, check_termination, check_infeasibility, scaling, MIN_SCALING, MAX_SCALING, adaptive_rho, adaptive_rho_interval, adaptive_rho_tolerance, verbose_timing, RHO_MIN, RHO_MAX, RHO_TOL, decompose, complete_dual, time_limit, obj_true, obj_true_tol, merge_strategy, colo_transformation)
end
end
