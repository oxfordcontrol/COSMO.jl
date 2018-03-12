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

  # constraint l <= Ax <= u
  mL = length(lInd)
  mU = length(uInd)

  # in this case also A has to be changed to get rid of -Inf and Inf values in
  Aa = [A[uInd,:] speye(mU) spzeros(mU,mL); A[lInd,:] spzeros(mL,mU) -speye(mL)]
  ba = [u[uInd];l[lInd]]
  nSlack = mL+mU
  K = OSSDPTypes.Cone(n,nSlack,[],[])
  if typeof(P) == Array{Float64,2}
    P = sparse(P)
    A = sparse(A)
  end
  Pa = blkdiag(P,spzeros(nSlack,nSlack))
  qa = [q;zeros(nSlack)]

  return Pa,qa,r, Aa,ba, K
  end


end # MODULE