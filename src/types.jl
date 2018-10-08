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

Object returned by the QOCS solver after calling `optimize!(model,settings)`. It has the following fields:


Fieldname | Type | Description
---  | --- | ---
x | Vector{T}| Primal variable
y | Vector{T}| Dual variable
s | Vector{T}| (Primal) set variable
objVal | T | Objective value
iter | Int64 | Number of iterations
status | Symbol | Solution status
info | QOCS.ResultInfo | Struct with more information
times | QOCS.ResultTimes | Struct with several measured times
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
    print(io,">>> QOCS - Results\nStatus: $(obj.status)\nIterations: $(obj.iter)\nOptimal Objective: $(round.(obj.objVal,digits=2))\nRuntime: $(round.(obj.times.solverTime*1000,digits=2))ms\nSetup Time: $(round.(obj.times.setupTime*1000,digits=2))ms\n")
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

Initializes an empty QOCS model that can be filled with problem data using `assemble!(model,P,q,constraints)`.
"""
mutable struct Model{T<:AbstractFloat}
    P::AbstractMatrix{T}
    q::Vector{T}
    A::AbstractMatrix{T}
    b::Vector{T}
    C::AbstractConvexSet{T}
    x0::Vector{T}
    y0::Vector{T}
    F #LDL Factorization
    m::Integer
    n::Integer
    flags::Flags

    function Model{T}() where{T}
        return new(
        spzeros(T,1,1),             #P
        T[],                        #q
        spzeros(T,1,1),             #A
        T[],                        #b
        ZeroSet{T}(1),              #C
        T[],                        #x0
        T[],                        #y0
        nothing,                    #F
        0,                          #m
        0,                          #n
        Flags())                    #Flags
    end
end

Model(args...) = Model{DefaultFloat}(args...)

# -------------------------------------
# Structure of internal iterate variables
# -------------------------------------

struct Variables{T}
    x::Vector{T}
    s::SplitVector{T}
    ν::Vector{T}
    μ::Vector{T}

    function Variables{T}(m::Int, n::Int, set::AbstractConvexSet{T}) where{T}
        x = zeros(T,n)
        s = zeros(T,m)
        μ = SplitVector(zeros(T,m),sets)
        ν = zeros(T,m)
        new(x,s,μ,ν)
    end
end

Variables(args...) = Variables{DefaultFloat}()


# -------------------------------------
# Top level container for all solver data
# -------------------------------------

mutable struct Workspace
    p::Model
    sm::ScaleMatrices
    vars::Variables
    ρ::Float64
    ρVec::Vector{Float64}
    Info::Info
    times::ResultTimes
    #constructor
    function Workspace(p::QOCS.Model,sm::ScaleMatrices)
        vars = Variables(Float64,p.m,p.n,p.C)
        ws = new(p,sm,vars,0.,Float64[],Info([0.]),ResultTimes())
        # hand over warm starting variables
        length(p.x0) == n && (ws.vars.x = p.x0)
        length(p.y0) == m && (ws.vars.μ = -p.y0)
        return ws
    end
end
