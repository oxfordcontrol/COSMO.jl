module Helper

export isNumericallyPosSemDef, isNumericallySymmetric, findNonSymmetricComponent, findDifferentElements, reCreateSparseMatrix,duplicateSparsityPattern



function isNumericallyPosSemDef(X,atol)
  X = X./2
  X = X+X'

  F = eigfact(X)
  if size(find(x-> x<-atol,F[:values]), 1) == 0
    return true
  else
    return false
  end
end

function isNumericallySymmetric(X,atol)
  n = size(X,2)
  for i = 1:n-1, j = i+1:n
    if abs(X[i,j] - X[j,i]) >= atol
      return false
    end
  end
    return true
end


function findNonSymmetricComponent(X)
  for i = 2:size(X,1), j = 1:(i-1)
    if abs(X[i,j] - X[j,i]) > 0.0
      return i,j,abs(X[i,j] - X[j,i])
    end
  end
end

function findDifferentElements(A,B)
  if size(A) != size(B)
    error("Matrices are not the same size")
  end
  m,n = size(A)
  diffEl = Array[]
  for iii = 1:m, jjj = 1:n
    if A[iii,jjj] != B[iii,jjj]
      push!(diffEl,[iii,jjj])
    end
  end
  return diffEl
end

function duplicateSparsityPattern(A)
  m,n = size(A)
  B = zeros(m,n)
  for iii = 1:m, jjj = 1:n
    if A[iii,jjj] != 0
      B[iii,jjj] = 1
    end
  end
  return B

end

function reCreateSparseMatrix(A)
  rowInd = A.rowval
  colPtr = A.colptr
  val = A.nzval

  #compute column indizes
  colInd = zeros(Int64,length(rowInd))
  cval = 1
  for iii=2:length(colPtr)
    currentPtr = colPtr[iii]
    prevPtr = colPtr[iii-1]
    colInd[prevPtr:currentPtr - 1] = cval
    cval+=1
  end

  # sort rowInd and vals
  p = sortperm(rowInd)
  return sparse(rowInd,colInd,val,size(A,1),size(A,2))

end
end #MODULE