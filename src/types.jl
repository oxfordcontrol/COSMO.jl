# -------------------------------------
# Results and related sub structures
# -------------------------------------

mutable struct ResultTimes{T<:AbstractFloat}
    solverTime::T
    setupTime::T
    graphTime::T
    factorTime::T
    iterTime::T
    projTime::T
    postTime::T
end

ResultTimes{T}() where{T} = ResultTimes{T}(0.,0.,0.,0.,0.,0.,0.)
ResultTimes(T::Type = DefaultFloat) = ResultTimes{T}()

struct ResultInfo{T<:AbstractFloat}
    rPrim::T
    rDual::T
end

ResultInfo(rp,rd) = ResultInfo{DefaultFloat}(rp,rd)

"""
Result{T <: AbstractFloat}

Object returned by the COSMO solver after calling `optimize!(model,settings)`. It has the following fields:


Fieldname | Type | Description
---  | --- | ---
x | Vector{T}| Primal variable
y | Vector{T}| Dual variable
s | Vector{T}| (Primal) set variable
objVal | T | Objective value
iter | Int64 | Number of iterations
status | Symbol | Solution status
info | COSMO.ResultInfo | Struct with more information
times | COSMO.ResultTimes | Struct with several measured times
"""
struct Result{T<:AbstractFloat}
    x::Vector{T}
    y::Vector{T}
    s::SplitVector{T}
    objVal::T
    iter::Int64
    status::Symbol
    info::ResultInfo
    times::ResultTimes

end

function Base.show(io::IO, obj::Result)
    print(io,">>> COSMO - Results\nStatus: $(obj.status)\nIterations: $(obj.iter)\nOptimal Objective: $(round.(obj.objVal,digits=2))\nRuntime: $(round.(obj.times.solverTime*1000,digits=2))ms\nSetup Time: $(round.(obj.times.setupTime*1000,digits=2))ms\n")
    obj.times.iterTime != NaN && print("Avg Iter Time: $(round.((obj.times.iterTime/obj.iter)*1000,digits=2))ms")
end

struct Info{T<:AbstractFloat}
    rho_updates::Vector{T}
end

# -------------------------------------
# Problem scaling
# -------------------------------------

struct ScaleMatrices{Tf <: AbstractFloat}
    D::Union{UniformScaling{Bool},Diagonal{Tf,Vector{Tf}} }
    Dinv::Union{UniformScaling{Bool},Diagonal{Tf,Vector{Tf}} }
    E::Union{UniformScaling{Bool},Diagonal{Tf,Vector{Tf}} }
    Einv::Union{UniformScaling{Bool},Diagonal{Tf,Vector{Tf}} }
    c::Base.RefValue{Tf}
    cinv::Base.RefValue{Tf}
end

ScaleMatrices(args...) = ScaleMatrices{DefaultFloat}(args...)

ScaleMatrices{T}() where {T} = ScaleMatrices(I,I,I,I,Base.RefValue{T}(one(T)),Base.RefValue{T}(one(T)))

function ScaleMatrices{T}(m,n) where{T}
    D    = Diagonal(ones(T,n))
    Dinv = Diagonal(ones(T,n))
    E    = Diagonal(ones(T,m))
    Einv = Diagonal(ones(T,m))
    c    = Base.RefValue{T}(one(T))
    cinv = Base.RefValue{T}(one(T))
    ScaleMatrices(D,Dinv,E,Einv,c,cinv)
end

# -------------------------------------
# a collection of flags
# -------------------------------------

mutable struct Flags
    FACTOR_LHS::Bool
    INFEASIBILITY_CHECKS::Bool
    REVERSE_SCALE_PROBLEM_DATA::Bool
    Flags() = new(true,true,true)
end


# -------------------------------------
# Problem data
# -------------------------------------

"""
Model()

Initializes an empty COSMO model that can be filled with problem data using `assemble!(model,P,q,constraints)`.
"""
mutable struct Model{T<:Real}
    P::AbstractMatrix{T}
    q::Vector{T}
    A::AbstractMatrix{T}
    b::Vector{T}
    C::CompositeConvexSet{T}
    x0::Vector{T}
    y0::Vector{T}
    model_size::Array{Integer,1}
    F::SuiteSparse.CHOLMOD.Factor{T}
    flags::Flags
    M::SparseMatrixCSC{T}

    function Model{T}() where{T}
        return new(
        spzeros(T,1,1),             #P
        T[],                        #q
        spzeros(T,1,1),             #A
        T[],                        #b
        COSMO.CompositeConvexSet([COSMO.ZeroSet{T}(1)]),     #C
        T[],                        #x0
        T[],                        #y0
        [0;0],                      #model size
        ldlt(sparse(1.0I,1,1)),     #F
        Flags(),                    #Flags
        spzeros(T,0,0))             #M
    end
end

Model(args...) = Model{DefaultFloat}(args...)

# -------------------------------------
# Structure of internal iterate variables
# -------------------------------------

struct Variables{T}
    x::Vector{T}
    s::SplitVector{T}
    μ::Vector{T}

    function Variables{T}(m::Int, n::Int, C::AbstractConvexSet{T}) where{T}
        m == C.dim || throw(DimensionMismatch("set dimension is not m"))
        x = zeros(T,n)
        s = SplitVector(zeros(T,m),C)
        μ = zeros(T,m)
        new(x,s,μ)
    end
end

Variables(args...) = Variables{DefaultFloat}(args...)


# -------------------------------------
# Top level container for all solver data
# -------------------------------------

mutable struct Workspace{T}
    p::Model{T}
    sm::ScaleMatrices{T}
    vars::Variables{T}
    ρ::T
    ρVec::Vector{T}
    Info::Info
    times::ResultTimes{Float64} #Always 64 bit regardless of data type?
    #constructor
    function Workspace{T}(p::COSMO.Model{T},sm::ScaleMatrices{T}) where{T}
        m, n  = p.model_size
        vars = Variables{T}(m,n,p.C)
        ws = new(p,sm,vars,zero(T),T[],Info([zero(T)]),ResultTimes())
        # hand over warm starting variables
        length(p.x0) == n && (ws.vars.x[:] = p.x0[:])
        length(p.y0) == m && (ws.vars.μ[:] = -p.y0[:])
        return ws
    end
end

Workspace(args...) = Workspace{DefaultFloat}(args...)
