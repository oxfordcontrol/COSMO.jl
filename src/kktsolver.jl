using QDLDL, COSMO


# -------------------------------------
# abstract KKT Solver.  Should implement:
# refactor! and solve! methods
# -------------------------------------
abstract type AbstractKKTSolver end

# -------------------------------------
# KKT Solver using QDLDL package
# -------------------------------------
"""
    QdldlKKTSolver()

Creates a solver for the KKT system that uses the LDL factorisation in QDLDL.
"""
struct QdldlKKTSolver <: AbstractKKTSolver

    m::Integer
    n::Integer
    ldlfact::QDLDL.QDLDLFactorisation

    function QdldlKKTSolver(m,n,ldlfact)
        positive_inertia(ldlfact) == n || error("Objective function is not convex.")
        new(m,n,ldlfact)
    end

end


function QdldlKKTSolver(P,A,sigma,rho)
    J = _make_Jacobian(P,A,sigma,rho)
    QdldlKKTSolver(size(A,1),size(A,2),qdldl(J))
end

_make_Jacobian(P,A,sigma::AbstractFloat,rho::AbstractFloat) = [P+sigma*I A'; A (-1/rho)*I]
_make_Jacobian(P,A,sigma::AbstractArray,rho::AbstractArray) = [P+Diagonal(sigma) A'; A Diagonal(-1 ./ rho)]
_make_Jacobian(P,A,sigma::AbstractFloat,rho::AbstractArray) = [P+sigma*I A'; A Diagonal(-1 ./ rho)]


function update_rho!(s::QdldlKKTSolver, rho)
    QDLDL.update_diagonal!(s.ldlfact,(s.n+1):(s.n+s.m),(-1 ./ rho))
end


function solve!(s::QdldlKKTSolver, x, b)
    x .= b
    QDLDL.solve!(s.ldlfact,x)
end
