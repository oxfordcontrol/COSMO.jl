# -------------------------------------
# Results and related sub structures
# -------------------------------------

mutable struct ResultTimes{T <: AbstractFloat}
	solver_time::T
	setup_time::T
	graph_time::T
	factor_time::T
	iter_time::T
	proj_time::T
	post_time::T
end

ResultTimes{T}() where{T} = ResultTimes{T}(0., 0., 0., 0., 0., 0., 0.)
ResultTimes(T::Type = DefaultFloat) = ResultTimes{T}()

struct ResultInfo{T <: AbstractFloat}
	r_prim::T
	r_dual::T
end

ResultInfo(rp, rd) = ResultInfo{DefaultFloat}(rp, rd)

"""
	Result{T <: AbstractFloat}

Object returned by the COSMO solver after calling `optimize!(model)`. It has the following fields:


Fieldname | Type | Description
---  | --- | ---
x | Vector{T}| Primal variable
y | Vector{T}| Dual variable
s | Vector{T}| (Primal) set variable
obj_val | T | Objective value
iter | Int64 | Number of iterations
status | Symbol | Solution status
info | COSMO.ResultInfo | Struct with more information
times | COSMO.ResultTimes | Struct with several measured times
"""
struct Result{T <: AbstractFloat}
    x::Vector{T}
    y::Vector{T}
    s::Vector{T}
    obj_val::T
    iter::Int64
    status::Symbol
    info::ResultInfo{T}
    times::ResultTimes{T}

    function Result{T}() where {T <: AbstractFloat}
      return new(zeros(T, 1), zeros(T, 1), zeros(T, 1), zero(T), 0, :Unsolved, ResultInfo{T}(0.,0.), ResultTimes{T}())
    end

    function Result{T}(x, y, s, obj_val, iter, status, info, times) where {T <: AbstractFloat}
      return new(x, y, s, obj_val, iter, status, info, times)
    end

end

function Base.show(io::IO, obj::Result)
	print(io,">>> COSMO - Results\nStatus: $(obj.status)\nIterations: $(obj.iter)\nOptimal Objective: $(round.(obj.obj_val, digits = 2))\nRuntime: $(round.(obj.times.solver_time * 1000, digits = 2))ms\nSetup Time: $(round.(obj.times.setup_time * 1000, digits = 2))ms\n")
	obj.times.iter_time != NaN && print("Avg Iter Time: $(round.((obj.times.iter_time / obj.iter) * 1000, digits = 2))ms")
end

struct Info{T <: AbstractFloat}
	rho_updates::Vector{T}
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

function ScaleMatrices{T}(m, n) where{T}
	D    = Diagonal(ones(T, n))
	Dinv = Diagonal(ones(T, n))
	E    = Diagonal(ones(T, m))
	Einv = Diagonal(ones(T, m))
	c    = Base.RefValue{T}(one(T))
	cinv = Base.RefValue{T}(one(T))
	ScaleMatrices(D, Dinv, E, Einv, c, cinv)
end

# -------------------------------------
# a collection of flags
# -------------------------------------

mutable struct Flags
	FACTOR_LHS::Bool
	INFEASIBILITY_CHECKS::Bool
	REVERSE_SCALE_PROBLEM_DATA::Bool
	Flags() = new(true, true, true)
end

# -------------------------------------
# Problem data
# -------------------------------------

mutable struct ProblemData{T<:Real}
	P::AbstractMatrix{T}
	q::Vector{T}
	A::AbstractMatrix{T}
	b::Vector{T}
	C::CompositeConvexSet{T}
	model_size::Array{Integer,1}

	function ProblemData{T}() where{T}
		return new(
			spzeros(T, 1, 1),             #P
			T[],                        #q
			spzeros(T, 1, 1),             #A
			T[],                        #b
			COSMO.CompositeConvexSet([COSMO.ZeroSet{T}(1)]),     #C
			[0; 0])                 #model size
	end
end

ProblemData(args...) = ProblemData{DefaultFloat}(args...)

# -------------------------------------
# Structure of internal iterate variables
# -------------------------------------

struct Variables{T}
	x::Vector{T}
	s::SplitVector{T}
	μ::Vector{T}

	function Variables{T}(m::Int, n::Int, C::AbstractConvexSet{T}) where{T}
		m == C.dim || throw(DimensionMismatch("set dimension is not m"))
		x = zeros(T, n)
		s = SplitVector(zeros(T, m), C)
		μ = zeros(T, m)
		new(x, s, μ)
	end
end

Variables(args...) = Variables{DefaultFloat}(args...)


# -------------------------------------
# Top level container for all solver data
# -------------------------------------
"""
	Workspace()

Initializes an empty COSMO model that can be filled with problem data using `assemble!(model, P, q,constraints, [settings])`.
"""
mutable struct Workspace{T}
	p::ProblemData{T}
	settings::Settings
	sm::ScaleMatrices{T}
	vars::Variables{T}
	ρ::T
	ρvec::Vector{T}
	F::SuiteSparse.CHOLMOD.Factor{T}
	M::SparseMatrixCSC{T}
	flags::Flags
	Info::Info
	times::ResultTimes{Float64} #Always 64 bit regardless of data type?
	#constructor
	function Workspace{T}() where {T}
		p = ProblemData{T}()
		sm = ScaleMatrices{T}()
		vars = Variables{T}(1, 1, p.C)
		return new(p, Settings(), sm, vars, zero(T), T[], ldlt(sparse(1.0I, 1, 1)), spzeros(0, 0), Flags(), Info([zero(T)]), ResultTimes())
	end
end
Workspace(args...) = Workspace{DefaultFloat}(args...)

# Type alias facing the user
"""
	Model()

Initializes an empty COSMO model that can be filled with problem data using `assemble!(model, P, q,constraints, [settings])`.
"""
const Model = Workspace;
