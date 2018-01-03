export create_scs_matrix, create_scs_settings, create_scs_data, create_scs_cone


# Takes a vector or matrix or sparse matrix A and creates an SCSMatrix
function create_scs_matrix(m::Int, n::Int, A::SCSVecOrMatOrSparse)
    A_sparse = sparse(A)

    values = A_sparse.nzval * 1.0
    rowval = convert(Array{Int, 1}, A_sparse.rowval .- 1)
    colptr = convert(Array{Int, 1}, A_sparse.colptr .- 1)

    return SCSMatrix(pointer(values), pointer(rowval), pointer(colptr), m, n)
end

function create_scs_settings(normalize=1::Int, scale=converter(Cdouble, 5.0)::Cdouble, rho_x=convert(Cdouble,1e-3)::Cdouble,
                        max_iters=2500::Int, eps=converter(Cdouble, 1e-3)::Cdouble, alpha=convert(Cdouble, 1.8)::Cdouble,
                        cg_rate=convert(Cdouble,2)::Cdouble, verbose=1::Int, warm_start=0::Int)
    return SCSSettings(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate, verbose, warm_start)
end

# Create an SCSData type
# We assume we are solving a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# A is the matrix with m rows and n cols
# b is of length m x 1
# c is of length n x 1
# refer to create_scs_cone for K
function create_scs_data(;m::Int=nothing, n::Int=nothing, A::Ptr{SCSMatrix}=nothing,
        b::Ptr{Cdouble}=nothing,  c::Ptr{Cdouble}=nothing, max_iters=20000::Int,
        eps=convert(Cdouble, 1e-4)::Cdouble, alpha=convert(Cdouble, 1.8)::Cdouble,
        rho_x=convert(Cdouble, 1e-3)::Cdouble, scale=convert(Cdouble, 5.0)::Cdouble,
        cg_rate=convert(Cdouble, 2)::Cdouble, verbose=1::Int,
        normalize=1::Int, warm_start=0::Int, options...)

    for (k, v) in options
        @eval(($k) = ($v))
    end
    stgs = create_scs_settings(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate, verbose, warm_start)
    return SCSData(m, n, A, b, c, pointer([stgs]))
end


# Refer to comment above
function create_scs_data(m::Int, n::Int, A::Ptr{SCSMatrix}, b::Ptr{Cdouble}, c::Ptr{Cdouble}; options...)
    return create_scs_data(m=m, n=n, A=A, b=b, c=c; options...)
end


# Refer to comment above
function create_scs_data(m::Int, n::Int, A::SCSVecOrMatOrSparse, b::Array{Float64,},
        c::Array{Float64,}; options...)
    if size(b, 1) != m || size(b, 2) != 1 || size(c, 1) != n || size(c, 2) != 1
        error("Size of b must be m x 1 and size of c must be n x 1")
    end
    A = [create_scs_matrix(m, n, A)]
    return create_scs_data(m=m, n=n, A=pointer(A), b=pointer(b), c=pointer(c); options...)
end


# Create an SCSCone type
# We assume we are solving a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
#
# f (num primal zero / dual free cones, i.e. primal equality constraints)
# l (num linear cones)
# q (array of SOCs sizes)
# s (array of SDCs sizes)
# ep (num primal exponential cones)
# ed (num dual exponential cones).
function create_scs_cone(f::Int, l::Int, q::Ptr{Int}, qsize::Int, s::Ptr{Int},
        ssize::Int, ep::Int, ed::Int, p::Ptr{Float64}, psize::Int)
    return SCSCone(f, l, q, qsize, s, ssize, ep, ed, p, psize)
end

# Refer to comment above
function create_scs_cone(f::Int, l::Int, q::Array{Int,}, qsize::Int, s::Array{Int,},
        ssize::Int, ep::Int, ed::Int, p::Array{Float64,}, psize::Int)
    return SCSCone(f, l, pointer(q), qsize, pointer(s), ssize, ep, ed, pointer(p), psize)
end


# We assume we are solving a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
#
#
# A is the matrix with m rows and n cols
# b is of length m x 1
# c is of length n x 1
#
# f (num primal zero / dual free cones, i.e. primal equality constraints)
# l (num linear cones)
# q (array of SOCs sizes)
# s (array of SDCs sizes)
# ep (num primal exponential cones)
# ed (num dual exponential cones).
#
# Returns object of type Solution
# type Solution with
# x, y, s, status (ASCII string), ret_val (numerical status)
#
function SCS_solve(m::Int, n::Int, A::SCSVecOrMatOrSparse, b::Array{Float64,},
        c::Array{Float64,}, f::Int, l::Int, q::Array{Int,}, qsize::Int, s::Array{Int,},
        ssize::Int, ep::Int, ed::Int, p::Array{Float64,}, psize::Int,
        primal_sol::Vector{Float64}=Float64[],
        dual_sol::Vector{Float64}=Float64[],
        slack::Vector{Float64}=Float64[];
        options...)

    data = create_scs_data(m, n, A, b, c; options...)
    cone = create_scs_cone(f, l, q, qsize, s, ssize, ep, ed, p, psize)

    if (:warm_start, true) in options && length(primal_sol) == n && length(dual_sol) == m && length(slack) == m
        x = primal_sol
        y = dual_sol
        s = slack
    else
        x = zeros(n)
        y = zeros(m)
        s = zeros(m)
    end
    solution = SCSSolution(pointer(x), pointer(y), pointer(s))
    status, solution, info, p_work = SCS_solve(data, cone, solution)
    SCS_finish(p_work)
    return Solution(x, y, s, status)

end

# for legacy API maintenance, sets power cones to zero
function SCS_solve(m::Int, n::Int, A::SCSVecOrMatOrSparse, b::Array{Float64,},
        c::Array{Float64,}, f::Int, l::Int, q::Array{Int,}, qsize::Int, s::Array{Int,},
        ssize::Int, ep::Int, ed::Int,
        primal_sol::Vector{Float64}=Float64[],
        dual_sol::Vector{Float64}=Float64[],
        slack::Vector{Float64}=Float64[];
        options...)
        return SCS_solve(m, n, A, b, c, f, l, q, qsize, s, ssize,
        ep, ed, Float64[], 0, primal_sol, dual_sol, slack; options...)
end