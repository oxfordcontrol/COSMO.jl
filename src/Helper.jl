module Helper

export isNumericallyPosSemDef, isNumericallySymmetric, findNonSymmetricComponent



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
      return i,j
    end
  end
end

end #MODULE