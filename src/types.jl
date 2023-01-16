# -------------------------------------
# Results and related sub structures
# -------------------------------------

"""
	ResultTimes

Part of the Result object returned by the solver. ResultTimes contains timing results for certain parts of the algorithm:

Time Name  | Description
---  | :---
solver_time | Total time used to solve the problem
setup_time |  Setup time = graph\\_time + init\\_factor\\_time + scaling\\_time
scaling_time | Time to scale the problem data
graph_time | Time used to perform chordal decomposition
init\\_factor\\_time | Time used for initial factorisation of the system of linear equations
factor\\_update\\_time | Sum of times used to refactor the system of linear equations due to rho
iter_time | Time spent in iteration loop
proj_time | Time spent in projection functions
post_time | Time used for post processing
update_time | Time spent in the update! function of the accelerator
accelerate_time | Time spent in the accelerate! function of the accelerator

By default COSMO only measures `solver_time`, `setup_time` and `proj_time`. To measure the other times set `verbose_timing = true`.
"""
mutable struct ResultTimes{T <: AbstractFloat}
	solver_time::T
	setup_time::T
	scaling_time::T
	graph_time::T
	init_factor_time::T
	factor_update_time::T
	iter_time::T
	proj_time::T
	post_time::T
	update_time::T
	accelerate_time::T
end

ResultTimes{T}() where{T} = ResultTimes{T}(T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN), T(NaN))
ResultTimes(T::Type = DefaultFloat) = ResultTimes{T}()

function Base.show(io::IO, obj::ResultTimes)
  obj.iter_time != 0 ? verbose = true : verbose = false
  print(io,"Solver time:\t$(round.(obj.solver_time, digits = 4))s ($(round.(obj.solver_time * 1000, digits = 2))ms)\n",
"Setup time:\t$(round.(obj.setup_time, digits = 4))s ($(round.(obj.setup_time * 1000, digits = 2))ms)\n",
"Proj time:\t$(round.(obj.proj_time, digits = 4))s ($(round.(obj.proj_time * 1000, digits = 2))ms)\n")
  if verbose
    print(io,"Iter time:\t$(round.(obj.iter_time, digits = 4))s ($(round.(obj.iter_time * 1000, digits = 2))ms)\n",
    "Scaling time:\t$(round.(obj.scaling_time, digits = 4))s ($(round.(obj.scaling_time * 1000, digits = 2))ms)\n",
    "Graph time:\t$(round.(obj.graph_time, digits = 4))s ($(round.(obj.graph_time * 1000, digits = 2))ms)\n",
    "Initial Factor time:\t$(round.(obj.init_factor_time, digits = 4))s ($(round.(obj.init_factor_time * 1000, digits = 2))ms)\n",
    "Factor update time:\t$(round.(obj.factor_update_time, digits = 4))s ($(round.(obj.factor_update_time * 1000, digits = 2))ms)\n",
    "Post time:\t$(round.(obj.post_time, digits = 4))s ($(round.(obj.post_time * 1000, digits = 2))ms)\n",
    "Accelerator - update time:\t$(round.(obj.update_time, digits = 4))s ($(round.(obj.update_time * 1000, digits = 2))ms)\n",
    "Accelerator - accelerate time:\t$(round.(obj.accelerate_time, digits = 4))s ($(round.(obj.accelerate_time * 1000, digits = 2))ms)\n")
  end
end

"""
    ResultInfo{T <: AbstractFloat}

Object that contains further information about the primal residual, the dual residuals and the rho updates.
"""
struct ResultInfo{T <: AbstractFloat}
	r_prim::T
	r_dual::T
	max_norm_prim::T
	max_norm_dual::T
	rho_updates::Vector{T}
end

ResultInfo(rp, rd, ro, rho_updates) = ResultInfo{DefaultFloat}(rp, rd, rho_updates)

"""
    Result{T <: AbstractFloat}

Object returned by the COSMO solver after calling `optimize!(model)`. It has the following fields:


Fieldname | Type | Description
---  | :--- | :---
x | Vector{T}| Primal variable
y | Vector{T}| Dual variable
s | Vector{T}| (Primal) set variable
obj_val | T | Objective value
iter | Int | Total number of ADMM iterations (incl. safeguarding_iter)
safeguarding_iter | Int | Number of iterations due to safeguarding of accelerator
status | Symbol | Solution status
info | COSMO.ResultInfo | Struct with more information
times | COSMO.ResultTimes | Struct with several measured times
"""
struct Result{T <: AbstractFloat}
    x::Vector{T}
    y::Vector{T}
    s::Vector{T}
    obj_val::T
    iter::Int
    safeguarding_iter::Int
    status::Symbol
    info::ResultInfo
    times::ResultTimes

    function Result{T}() where {T <: AbstractFloat}
      return new(zeros(T, 1), zeros(T, 1), zeros(T, 1), zero(T), 0, 0, :Unsolved, ResultInfo{T}(0., 0., 0., 0., T[]), ResultTimes{T}())
    end

    function Result{T}(x, y, s, obj_val, iter, safeguarding_iter, status, info, times) where {T <: AbstractFloat}
      return new(x, y, s, obj_val, iter, safeguarding_iter, status, info, times)
    end

end

function Base.show(io::IO, obj::Result)
	print(io,">>> COSMO - Results\nStatus: ")
	if obj.status == :Solved
		result_color = :green
	else
		result_color = :red
	end
	printstyled("$(obj.status)\n", color = result_color)
	println("Iterations: $(obj.iter) (incl. $(obj.safeguarding_iter) safeguarding iterations)\nOptimal Objective: $(@sprintf("%.4g", obj.obj_val))\nRuntime: $(round.(obj.times.solver_time * 1000, digits = 2))ms\nSetup Time: $(round.(obj.times.setup_time * 1000, digits = 2))ms\n")
	!isnan(obj.times.iter_time) && print("Avg Iter Time: $(round.((obj.times.iter_time / obj.iter) * 1000, digits = 2))ms")
end

# -------------------------------------
# Problem scaling
# -------------------------------------

struct ScaleMatrices{Tf <: AbstractFloat}
	D::Union{UniformScaling{Bool}, Diagonal{Tf, Vector{Tf}} }
	Dinv::Union{UniformScaling{Bool}, Diagonal{Tf, Vector{Tf}} }
	E::Union{UniformScaling{Bool}, Diagonal{Tf, Vector{Tf}} }
	Einv::Union{UniformScaling{Bool}, Diagonal{Tf, Vector{Tf}} }
	c::Base.RefValue{Tf}
	cinv::Base.RefValue{Tf}
end

ScaleMatrices(args...) = ScaleMatrices{DefaultFloat}(args...)

ScaleMatrices{T}() where {T} = ScaleMatrices(I, I, I, I, Base.RefValue{T}(one(T)), Base.RefValue{T}(one(T)))

function ScaleMatrices{T}(m, n) where {T <: AbstractFloat}
	D    = Diagonal(ones(T, n))
	Dinv = Diagonal(ones(T, n))
	E    = Diagonal(ones(T, m))
	Einv = Diagonal(ones(T, m))
	c    = Base.RefValue{T}(one(T))
	cinv = Base.RefValue{T}(one(T))
	ScaleMatrices(D, Dinv, E, Einv, c, cinv)
end


# -------------------------------------
# Problem data
# -------------------------------------

mutable struct ProblemData{T <: AbstractFloat}
	P::AbstractMatrix{T}
	q::Vector{T}
	A::AbstractMatrix{T}
	b::Vector{T}
	C::CompositeConvexSet{T}
	model_size::Array{Integer,1}

	function ProblemData{T}() where {T <: AbstractFloat}
		return new(
			spzeros(T, 1, 1),             #P
			T[],                        #q
			spzeros(T, 1, 1),             #A
			T[],                        #b
			COSMO.CompositeConvexSet{T}([COSMO.ZeroSet{T}(1)]),     #C
			[0; 0])                 #model size
	end
end

ProblemData(args...) = ProblemData{DefaultFloat}(args...)

# ---------------------------
# Struct to hold clique and sparsity data for a constraint
# ---------------------------

mutable struct SparsityPattern
  sntree::SuperNodeTree
  ordering::Array{Int}
  reverse_ordering::Array{Int}
  row_range::UnitRange{Int} # the starting row of the psd cone in the original problem
  cone_ind::Int # this is the ind of the original psd cone in ws.p.C that is decomposed
  nz_ind_map::SparseVector{Int, Int} # maps a matrix entry k = svec(i, j) to the location of the entry in the sparse data structure

  # constructor for sparsity pattern
  function SparsityPattern(L::SparseMatrixCSC, N::Int, ordering::Array{Int, 1}, merge_strategy::Union{Type{<: AbstractMergeStrategy}, OptionsFactory{<: AbstractMergeStrategy}}, row_range::UnitRange{Int}, cone_ind::Int, nz_ind_map::SparseVector{Int, Int})

    merge_strategy = merge_strategy()
    sntree = SuperNodeTree(L, merge_strategy)

    # clique merging only if more than one clique present
    if sntree.num > 1
		merge_cliques!(sntree)
	elseif merge_strategy isa AbstractGraphBasedMerge
		# if no merging attempt happens for a clique graph, we still have to convert the snd and sep back to Array{Array{Int}, 1}
		# for consistency
		sntree.snd = sort.(collect.(sntree.snd))
		sntree.sep = sort.(collect.(sntree.sep))
	end
    # reorder vertices in supernodes to have consecutive order
    # necessary for equal column structure for psd completion
    reorder_snd_consecutively!(sntree, ordering)


    # for each clique determine the number of entries of the block represented by that clique
    calculate_block_dimensions!(sntree)

    return new(sntree, ordering, invperm(ordering), row_range, cone_ind, nz_ind_map)
  end

  # For debugging
  function SparsityPattern(sntree::SuperNodeTree, ordering::Array{Int}, reverse_ordering::Array{Int}, row_range::UnitRange{Int}, cone_ind::Int)
    return new(sntree, ordering, reverse_ordering, row_range, cone_ind)
  end
end

# -------------------------------------
# Chordal Decomposition Information
# -------------------------------------
mutable struct ChordalInfo{T <: AbstractFloat}
  decompose::Bool # an internal flag to check if problem has been decomposed
  originalM::Int
  originalN::Int
  originalC::CompositeConvexSet{T}
  H::SparseMatrixCSC{T}
  sp_arr::Array{COSMO.SparsityPattern}
  psd_cones_ind::Array{Int} # stores the position of decomposable psd cones in the composite convex set
  num_psd_cones::Int # number of psd cones of original problem
  num_decomposable::Int #number of decomposable cones
  num_decom_psd_cones::Int #total number of psd cones after decomposition
  L::SparseMatrixCSC{T} #pre allocate memory for QDLDL
  cone_map::Dict{Int, Int} # map every cone in the decomposed problem to the equivalent or undecomposed cone in the original problem

  function ChordalInfo{T}(problem::COSMO.ProblemData{T}, settings::COSMO.Settings) where {T <: AbstractFloat}
    originalM = problem.model_size[1]
    originalN = problem.model_size[2]
    originalC = deepcopy(problem.C)
    num_psd_cones = length(findall(x -> typeof(x) <: Union{PsdConeTriangle{T}, PsdCone{T}} , problem.C.sets))
    # allocate sparsity pattern for each cone
    sp_arr = Array{COSMO.SparsityPattern}(undef, num_psd_cones)
    cone_map = Dict{Int, Int}()

    return new(settings.decompose, originalM, originalN, originalC, spzeros(1, 1), sp_arr, Int[], num_psd_cones, 0, 0, spzeros(1, 1), cone_map)
  end

	function ChordalInfo{T}() where{T}
		C = COSMO.CompositeConvexSet{T}([COSMO.ZeroSet{T}(1)])
		return new(false, 0, 0, C, spzeros(1, 1), COSMO.SparsityPattern[], [1])
	end

end

# -------------------------------------
# Structure of internal iterate variables
# -------------------------------------

struct Variables{T}
	w::Vector{T}
	w_prev::Vector{T}
	x::SubArray{T}
	s::SplitVector{T}
	μ::Vector{T}

	function Variables{T}(m::Int, n::Int, C::AbstractConvexSet{T}) where {T <: AbstractFloat}
		m == C.dim || throw(DimensionMismatch("set dimension is not m"))
		w = zeros(T, n + m)
		w_prev = zeros(T, n + m)
		x = view(w_prev, 1:n) #x is a view onto w_prev, to keep it in sync with s and μ 
		s = SplitVector(zeros(T, m), C)
		μ = zeros(T, m)
		new(w, w_prev, x, s, μ)
	end
end

Variables(args...) = Variables{DefaultFloat}(args...)

mutable struct IterateHistory
	x_data::AbstractMatrix
	s_data::AbstractMatrix
	y_data::AbstractMatrix
	v_data::AbstractMatrix
	r_prim_data::AbstractVector
	r_dual_data::AbstractVector
	eta_data::AbstractArray
	alpha_data::AbstractArray
	cond_data::AbstractVector
	aa_fail_data::AbstractVector

	function IterateHistory(m, n, mem)
		new(zeros(n, 0), zeros(m, 0), zeros(m, 0), zeros(m + n, 0), Float64[], Float64[], zeros(mem, 0), zeros(mem + 1, 0), Float64[], Int64[])
	end
end

function update_iterate_history!(history::IterateHistory, x, s, y, v, r_prim, r_dual, eta::Vector{Float64}, cond::Float64)
	history.x_data = hcat(history.x_data, x)
	history.s_data = hcat(history.s_data, s)
	history.y_data = hcat(history.y_data, y)
	history.v_data = hcat(history.v_data, v)
	history.cond_data = push!(history.cond_data, cond)
	push!(history.r_prim_data, r_prim)
	push!(history.r_dual_data, r_dual)

	# alphas = compute_alphas(eta)
	# history.alpha_data = hcat(history.alpha_data, alphas)
	# history.eta_data = hcat(history.eta_data, eta)
end

struct UtilityVariables{T}
  vec_m::Vector{T}
  vec_n::Vector{T}
  vec_n2::Vector{T}

  function UtilityVariables{T}(m::Int, n::Int) where {T <: AbstractFloat}
    new(zeros(T, m), zeros(T, n), zeros(T, n))
  end
end

UtilityVariables(args...) = UtilityVariables{DefaultFloat}(args...)

# -------------------------------------
# a collection of state flags
# -------------------------------------

mutable struct States
	IS_ASSEMBLED::Bool # the workspace has been assembled with problem data
	IS_OPTIMIZED::Bool # the optimization function has been called on the model
	IS_CHORDAL_DECOMPOSED::Bool # the problem has been decomposed
	KKT_FACTORED::Bool # the KKT matrix has been factored
	IS_SCALED::Bool # the problem data has been scaled
	States() = new(false, false, false, false, false)
end


# -------------------------------------
# Top level container for all solver data
# -------------------------------------
"""
	Workspace{T <: AbstractFloat}()

Initializes an empty COSMO model that can be filled with problem data using `assemble!(model, P, q,constraints; [settings, x0, s0, y0])`.
"""
mutable struct Workspace{T}
	p::ProblemData{T}
	settings::Settings{T}
	sm::ScaleMatrices{T}
	ci::ChordalInfo{T}
	vars::Variables{T}
  	utility_vars::UtilityVariables{T}
	δx::Vector{T}
	δy::SplitVector{T}
	s_tl::Vector{T}
	ls::Vector{T}
	sol::Vector{T}
	x_tl::SubArray{T}
	ν::SubArray{T}
	ρ::T
	ρvec::Vector{T}
	kkt_solver::Union{AbstractKKTSolver,Nothing}
	states::States
	rho_updates::Vector{T} #keep track of the rho updates and the number of refactorisations
	rho_update_due::Bool # a flag that indicates that at the next possible iteration ρ should be updated
	infeasibility_check_due::Bool # a flag that indicates that at the next possible iteration we should check for infeasibility
	times::ResultTimes{Float64} #always 64 bit regardless of data type
	row_ranges::Array{UnitRange{Int}, 1} # store a set_ind -> row_range map
	safeguarding_iter::Int # count extra ADMM iterations due to bad quality accelerated steps
	accelerator::AbstractAccelerator
	accelerator_active::Bool
	activation_reason::AbstractActivationReason
	#constructor
	function Workspace{T}() where {T <: AbstractFloat}
		p = ProblemData{T}()
		sm = ScaleMatrices{T}()
		vars = Variables{T}(1, 1, p.C)
    	uvars = UtilityVariables{T}(1, 1)
		ci = ChordalInfo{T}()
		δx = zeros(T, 0)
		δy = SplitVector(zeros(T, 1), p.C)
		s_tl = zeros(T,0)
		ls = zeros(T, 0)
		sol = zeros(T, 1)
		x_tl = view(sol, 1:1)
		ν = view(sol, 1:1)
		return new(p, Settings{T}(), sm, ci, vars,  uvars, δx, δy, s_tl, ls, sol, x_tl, ν, zero(T), T[], nothing, States(), T[], false, false, ResultTimes(), [0:0], 0, EmptyAccelerator(), false, ImmediateActivation())
	end
end
Workspace(args...) = Workspace{DefaultFloat}(args...)

Base.show(io::IO, model::COSMO.Workspace{T}) where {T} = println(io, "A COSMO Model with Float precision: $(T)")


# Type alias facing the user
"""
	Model{T <: AbstractFloat}()

Initializes an empty COSMO model that can be filled with problem data using `assemble!(model, P, q,constraints; [settings, x0, s0, y0])`.
"""
const Model = Workspace;
