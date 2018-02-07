module Helper

export isNumericallyPosSemDef, isNumericallySymmetric



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




end #MODULE