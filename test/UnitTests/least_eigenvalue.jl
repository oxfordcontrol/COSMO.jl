using COSMO, Test, LinearAlgebra, SparseArrays
# this problem computes the least eigenvalue of a Hermitian matrix c

# min dot(c, X)
# s.t.  tr(X) = 1, X ⪰ 0

tri(d) = div(d*(d+1), 2)

function mineig_raw(c::AbstractMatrix{R}) where {R}
    d = size(c, 1)
    T = real(R)
    vec_dim = R <: Complex ? d^2 : tri(d)    

    model = COSMO.Model{T}()
    vec_c = zeros(T, vec_dim)
    COSMO.extract_upper_triangle!(c, vec_c, sqrt(T(2)))

    P = zeros(T, vec_dim, vec_dim)

    id_vec = zeros(T, vec_dim)
    diagonal_indices = tri.(1:d)
    id_vec[diagonal_indices] .= 1

    cs1 = COSMO.Constraint(id_vec', -T(1), COSMO.ZeroSet)
    cs2 = COSMO.Constraint(Matrix(T(1)*I(vec_dim)), zeros(T,vec_dim), COSMO.PsdConeTriangle{T, R}(vec_dim))
    constraints = [cs1; cs2]

    assemble!(model, P, vec_c, constraints, settings = COSMO.Settings{T}(verbose = false, eps_abs = 1e-5, eps_rel = 1e-5))
    result = COSMO.optimize!(model)
    return result.obj_val
end

@testset "Complex PSD Cone" begin
    X = Hermitian(ComplexF64.([ 1  im 0;
                               -im 1 im;
                                0 -im 1]))
    @test mineig_raw(X) ≈ 1-sqrt(2) atol = 1e-4 rtol = 1e-4
end
