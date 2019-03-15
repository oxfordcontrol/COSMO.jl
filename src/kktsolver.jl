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

end


function QdldlKKTSolver(P,A,sigma::AbstractFloat,rho::AbstractFloat)
    J = [P+sigma*I A'; A (-1 ./ rho)*I]
    QdldlKKTSolver(size(A,1),size(A,2),qdldl(J))
end


function QdldlKKTSolver(P,A,sigma::Vector,rho::Vector)
    J = [P+Diagonal(sigma) A'; A Diagonal(-1 ./ rho)]
    QdldlKKTSolver(size(A,1),size(A,2),qdldl(J))
end


function update_rho!(s::QdldlKKTSolver, rho)
    QDLDL.update_diagonal!(s.ldlfact,(s.n+1):(s.n+s.m),(-1 ./ rho))
end


function solve!(s::QdldlKKTSolver, rhs)
    QDLDL.solve!(s.ldlfact,rhs)
end
