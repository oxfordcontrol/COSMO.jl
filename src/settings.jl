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

"Absctract supertype for accelerator activation."
abstract type AbstractActivationReason end

"Activate accelerator immediately."
struct ImmediateActivation <: AbstractActivationReason end


"""
	COSMO.Settings{T}(; kwargs) where {T <: AbstractFloat}

Creates a COSMO settings object that is used to pass user settings to the solver.

Argument | Description | Values (default)
--- | :--- | :---
rho | ADMM rho step | 0.1
sigma | ADMM sigma step | 1e-6
alpha | Relaxation parameter | 1.6
eps_abs | Absolute residual tolerance | 1e-5
eps_rel | Relative residual tolerance | 1e-5
eps\\_prim\\_inf | Primal infeasibility tolerance | 1e-5
eps\\_dual\\_inf | Dual infeasibility tolerance | 1e-5
max_iter | Maximum number of iterations | 5000
verbose | Verbose printing | false
verbose_timing | Verbose timing | false
kkt_solver | Linear System solver | `QdldlKKTSolver`
check_termination | Check termination interval | 25
check_infeasibility | Check infeasibility interval | 40
scaling | Number of scaling iterations | 10
adaptive_rho | Automatic adaptation of step size parameter | true
adaptive_rho_max_adaptions | Max number of rho adaptions | typemax(Int64) (deactivated)
decompose | Activate to decompose chordal psd constraints | true
complete_dual | Activate to complete the dual variable after decomposition | false
merge_strategy | Choose a strategy for clique merging | `CliqueGraphMerge`
compact_transformation | Choose how a decomposed problem is transformed | true
time_limit | Set solver time limit in s | 0 (deactivated)
accelerator | Acceleration scheme | `AndersonAccelerator{Type2}`
accelerator_activation | Accelerator activation | `ImmediateActivation`
safeguarding | Accelerator safeguarding | true
safeguarding_tol | Safeguarding tolerance | 2.0

"""
mutable struct Settings{T <: AbstractFloat}
	rho::T
	sigma::T
	alpha::T
	eps_abs::T
	eps_rel::T
	eps_prim_inf::T
	eps_dual_inf::T
	max_iter::Int
	verbose::Bool
	kkt_solver::Union{Type{<: AbstractKKTSolver}, OptionsFactory{<: AbstractKKTSolver}}
	check_termination::Int
	check_infeasibility::Int
	scaling::Int
	MIN_SCALING::T
	MAX_SCALING::T
	adaptive_rho::Bool
	adaptive_rho_interval::Int
	adaptive_rho_tolerance::T
	adaptive_rho_fraction::T
	adaptive_rho_max_adaptions::Int64
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
	accelerator::Union{Type{<: AbstractAccelerator}, OptionsFactory{<: AbstractAccelerator}}
	safeguarding::Bool
	safeguarding_tol::T
	#constructor
	function Settings{T}(;
		rho::Real=T(0.1),
		sigma::Real=T(1e-6),
		alpha::Real=T(1.6),
		eps_abs::Real=T(1e-5),
		eps_rel::Real=T(1e-5),
		eps_prim_inf::Real=T(1e-4),
		eps_dual_inf::Real=T(1e-4),
		max_iter::Integer=5000,
		verbose::Bool=false,
		kkt_solver=QdldlKKTSolver,
		check_termination::Int=25,
		check_infeasibility::Int=40,
		scaling::Integer=10,
		MIN_SCALING::Real = T(1e-4),
		MAX_SCALING::Real = T(1e4),
		adaptive_rho::Bool = true,
		adaptive_rho_interval::Int = 40,
		adaptive_rho_tolerance::Int = 5,
		adaptive_rho_fraction::Real = T(0.4),
		adaptive_rho_max_adaptions::Int = typemax(Int),
		verbose_timing::Bool = false,
		RHO_MIN::Real = T(1e-6),
		RHO_MAX::Real = T(1e6),
		RHO_TOL::Real = T(1e-4),
		RHO_EQ_OVER_RHO_INEQ::Real = T(1e3),
		COSMO_INFTY::Real = T(1e20),
		decompose::Bool = true,
    	complete_dual::Bool = false,
		time_limit::Real = zero(T),
		obj_true::Real = T(NaN),
		obj_true_tol::Real = T(1e-3),
		merge_strategy = CliqueGraphMerge,
		compact_transformation::Bool = true,
		accelerator = with_options(AndersonAccelerator{T, Type2{QRDecomp}, RestartedMemory, NoRegularizer}, mem = 10),
		safeguarding::Bool = true, 
		safeguarding_tol::T = T(2)
		) where {T <: AbstractFloat}
		if !isa(kkt_solver, OptionsFactory)
			kkt_solver = with_options(kkt_solver)
		end

		if !isa(merge_strategy, OptionsFactory)
			merge_strategy = with_options(merge_strategy)
		end

		if !isa(accelerator, OptionsFactory)
			accelerator = with_options(accelerator)
		end
		
		new(rho, sigma, alpha, eps_abs, eps_rel, eps_prim_inf, eps_dual_inf, max_iter, verbose, kkt_solver, check_termination, check_infeasibility, scaling, MIN_SCALING, MAX_SCALING, adaptive_rho, adaptive_rho_interval, adaptive_rho_tolerance, adaptive_rho_fraction, adaptive_rho_max_adaptions, verbose_timing, RHO_MIN, RHO_MAX, RHO_TOL, RHO_EQ_OVER_RHO_INEQ, COSMO_INFTY, decompose, complete_dual, time_limit, obj_true, obj_true_tol, merge_strategy, compact_transformation, accelerator, safeguarding, safeguarding_tol)
	end
end

# The default case it to return a Settings{Float64} object
Settings(args...; kwargs...) = Settings{DefaultFloat}(args...; kwargs...)


function Base.show(io::IO, obj::COSMO.Settings{T}) where {T <: AbstractFloat}
	println(io,"A COSMO.Settings{$(T)} object. To list the available options type `?` and `help?>COSMO.Settings`.")
end
eltype(::Type{<:Settings{T}}) where {T} = T


function Settings(d::Dict)

	# converts strings to julia types by look-up
	string_settings_converter = Dict("QdldlKKTSolver" => COSMO.QdldlKKTSolver, "CholmodKKTSolver" => COSMO.CholmodKKTSolver, "NoMerge" => COSMO.NoMerge, "ParentChildMerge" => COSMO.ParentChildMerge, "CliqueGraphMerge" => COSMO.CliqueGraphMerge)

	settings = COSMO.Settings{DefaultFloat}()
	# convert pure dictionary settings to COSMO.Settings object
	settings = COSMO.Settings()
	for (key, val) in d
		if key in ["kkt_solver"; "merge_strategy"]
			setfield!(settings, Symbol(key), with_options(string_settings_converter[val]))
		else
			setfield!(settings, Symbol(key), val)
		end
	end
	return settings
end
