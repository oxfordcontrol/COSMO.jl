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

  m = size(A,1)
  n = size(P,1)

  # Rewrite problem to OSSDP compatible format:
  # determine the indizes and sizes of box constraints with Inf or -Inf values
  uInd = find(x->(x < 1e19 ),u)
  lInd = find(x->(x > -1e19),l)

  # in this case also A has to be changed to get rid of -Inf and Inf values in
  Aa = [A[uInd,:]; -A[lInd,:]]
  ba = [u[uInd]; -l[lInd]]
  K = OSSDPTypes.Cone(0, size(Aa)[1], [], [])
  if typeof(Aa) == Array{Float64,2}
    Aa = sparse(Aa)
  end
  if typeof(P) == Array{Float64,2}
    P = sparse(P)
  end

  return P,q,r,Aa,ba,K
  end


end # MODULE