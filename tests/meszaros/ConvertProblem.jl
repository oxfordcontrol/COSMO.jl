# Module to convert the problem definitions of the Maros Meszaros' QP Test set into a OSSDP compatible format
# Main problem is to identify the rows that contain inequalities with Inf and -Inf, that have to be removed

module Converter

using OSSDPTypes
export convertProblem


function convertProblem(data)
  Q = data["Q"]
  A = data["A"]
  c = data["c"]
  ru = data["ru"]
  rl = data["rl"]
  lb = data["lb"]
  ub = data["ub"]

  m = size(A,1)
  n = size(Q,1)

  # Rewrite problem to OSSDP compatible format:
  P = Q
  q = c
  Aa = ba = K =  0

  # determine the indizes and sizes of box constraints with Inf or -Inf values
  ubInd = find(x->(x != Inf ),ub)
  lbInd = find(x->(x != -Inf),lb)
  rlInd = find(x->(x != -Inf),rl)
  ruInd = find(x->(x != Inf),ru)

  # constraint lb <= x <= ub
  In = eye(n)
  Iub = In[ubInd,:]
  Ilb = In[lbInd,:]
  mUB, = size(Iub)
  mLB, = size(Ilb)

  # constraint rl <= Ax <= ru
  Im = eye(m)
  Irl = Im[rlInd,:]
  Iru = Im[ruInd,:]
  mRL, = size(Irl)
  mRU, = size(Iru)

  # remove rows of A relating to Inf inequality constraints
  Au = A[ruInd,:]
  Al = A[rlInd,:]
  typ = 0
  # check if case Ax = b
  if norm(rl-ru,Inf) < 1e-4
    typ = 1
    b = rl
    Aa = [A zeros(m,mUB+mLB); Iub eye(mUB) zeros(mUB,mLB);Ilb zeros(mLB,mUB) -eye(mLB)]
    ba = [b;ub[ubInd];lb[lbInd]]
    nSlack = mUB+mLB
    K = OSSDPTypes.Cone(n,nSlack,[],[])
    Pa = blkdiag(P,spzeros(nSlack,nSlack))
    qa = [q;zeros(nSlack)]
  else
    typ = 2
    # in this case also A has to be changed to get rid of -Inf and Inf values in
    Aa = [Au eye(mRU) zeros(mRU,mRL+mUB+mLB); Al zeros(mRL,mRU) -eye(mRL) zeros(mRL,mLB+mUB); Iub zeros(mUB,mRL+mRU) eye(mUB) zeros(mUB,mLB);Ilb zeros(mLB,mRL+mRU+mUB) -eye(mLB)]
    ba = [ru[ruInd];rl[rlInd];ub[ubInd];lb[lbInd]]
    nSlack = mRL+mRU+mLB+mUB
    K = OSSDPTypes.Cone(n,nSlack,[],[])
    Pa = blkdiag(P,spzeros(nSlack,nSlack))
    qa = [q;zeros(nSlack)]
  end
  return Pa,qa,Aa,ba, K,typ
  end


end # MODULE