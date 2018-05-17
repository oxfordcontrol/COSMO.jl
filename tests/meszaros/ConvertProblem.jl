# Module to convert the problem definitions of the Maros Meszaros' QP Test set into a OSSDP compatible format
# Main problem is to identify the rows that contain inequalities with Inf and -Inf, that have to be removed

module Converter

using OSSDPTypes
export convertProblem


function convertProblem(data)
  P = data["P"]
  A = data["A"]
  q = data["q"]
  u = data["u"]
  l = data["l"]
  r = data["r"]

  # The matrix A for some reason segfaults in the QAFIRO problem.
  # We avoid this by recreating the sparse matrix
  I, J, V = findnz(A)
  A = sparse(I, J, V)

  # Rewrite problem to OSSDP compatible format:
  # determine the indizes and sizes of box constraints with Inf or -Inf values
  neq = find(u .!== l)
  eq = find(u .== l)
  Aeq = A[eq, :]
  beq = u[eq]

  A = A[neq, :]
  u = u[neq]
  l = l[neq]
  uInd = find(x->(x < 1e19 ),u)
  lInd = find(x->(x > -1e19),l)

  # in this case also A has to be changed to get rid of -Inf and Inf values in
  Aa = [Aeq; A[uInd,:]; -A[lInd,:]]
  ba = [beq; u[uInd]; -l[lInd]]
  K = OSSDPTypes.Cone(size(eq)[1], size(uInd)[1] + size(lInd)[1], [], [])
  if typeof(Aa) == Array{Float64,2}
    Aa = sparse(Aa)
  end
  if typeof(P) == Array{Float64,2}
    P = sparse(P)
  end

  return P,q,r,Aa,ba,K
  end


end # MODULE