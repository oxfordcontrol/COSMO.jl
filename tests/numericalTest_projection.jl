workspace()

rng = MersenneTwister(1234)

# create random symmetric matrix
X = rand(rng,5,5)
X = full(Symmetric(X))

# compute eigenvalue decomposition
F = eigfact(X)

Λ = diagm(F[:values])
Q = F[:vectors]

# set negative eigenvalues to 0
Xp = Q*max.(Λ,0)*Q'

println("Xp is posdef? $(isposdef(Xp))")

# compute eigenvalues
E = eigfact(Xp)
println("Eigenvalues: $(E[:values])")

# take eigenvector with negative eigenvalue
v = F[:vectors][:,2]

# compare the following results according to Paul
println("v'Xpv= $(v'*Xp*v)")
println("v'Q max(Λ,0)(Q'v) = $(v'*Q*max.(Λ,0)*Q'*v)")
