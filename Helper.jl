module Helper

export isNumericallyPosSemDef



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




end #MODULE