"""
Settings(;arg=val)

Creates a COSMO settings object that is used to pass user settings to the solver.

Argument | Description | Values (default)
--- | --- | ---
rho | ADMM rho step | 0.1
sigma | ADMM sigma step | 1e-6.
alpha | Relaxation parameter | 1.6
eps_abs | Absolute residual tolerance | 1e-4
eps_rel | Relative residual tolerance | 1e-4
eps_prim_inf | Primal infeasibility tolerance | 1e-4
eps_dual_inf | Dual infeasibility tolerance | 1e-4
max_iter | Maximum number of iterations | 2500
verbose | Verbose printing | false
verbose_timing | Verbose timing | false
check_termination | Check termination interval | 40
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
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
	time_limit::Float64
	obj_true::Float64
	obj_true_tol::Float64
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
		time_limit = 0.0,
		obj_true = NaN,
		obj_true_tol = 1e-3
		)
	new(rho, sigma, alpha, eps_abs, eps_rel, eps_prim_inf, eps_dual_inf, max_iter, verbose,  check_termination, check_infeasibility, scaling, MIN_SCALING, MAX_SCALING, adaptive_rho, adaptive_rho_interval, adaptive_rho_tolerance, verbose_timing, RHO_MIN, RHO_MAX, RHO_TOL, time_limit, obj_true, obj_true_tol)
end
end
