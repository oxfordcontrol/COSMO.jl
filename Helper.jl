module Helper

export isNumericallyPosSemDef



function isNumericallyPosSemDef(X,eps)
  X = X./2
  X = X+X'

  F = eigfact(X)
  if size(find(x-> x<eps,F[:values]), 1) == 0
    return true
  else
    return false
  end
end




end #MODULE