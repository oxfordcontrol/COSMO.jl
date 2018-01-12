# Test file to determine the influence of the relaxation parameter
# alpha onto the convergence of the algorithm

using PyPlot

# Define again the example problem
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0
# True solution obtained with ConvexJL + SCS
XTrue =[1.0559 0.369196 0.86829; 0.369196 0.129089 0.303588; 0.868288 0.303592 0.714004]
objTrue = 13.902243096458092


# Reformulate data to suit OSSDP input
P = zeros(3,3)
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]

αArr = 0.5:0.1:1.6
#αArr = 1.0:0.1:1.1
nn = 2000
iterToConv = zeros(length(αArr))
# test algorithm's convergence for different relaxation values
counter = 1
for α in αArr

  settings = sdpSettings(rho=1.0,sigma=1.0,alpha=α,max_iter=nn,verbose=false)

  # solve SDP problem
  res,dbg = solveSDP(P,q,A,b,settings)

  # Check after how many iterations a reasonable convergence to the true solution is achieved
  for i=1:nn
    X = reshape(dbg.x[i,:],3,3)
    S = reshape(dbg.x[i,:],3,3)
    cost = dbg.cost[i]
    F = eigfact(X)
    XλMin = minimum(F[:values])
    if norm(X-XTrue,Inf) < 1e-3 && norm(cost-objTrue)< 1e-3 &&  norm(vec(X)-vec(S),Inf)<1e-3 && XλMin > -1e-9
      iterToConv[counter] = i
      break
    else
      iterToConv[counter] = i
    end
  end
counter += 1
end

#plot debugging variables

PyPlot.figure(1,facecolor="white")
PyPlot.plot(αArr,iterToConv, linewidth=2.0)
PyPlot.grid(true)
PyPlot.xticks(αArr)
PyPlot.xlabel(L"α")
PyPlot.ylabel("Iterations to Convergence")
PyPlot.title("Influence of relaxation parameter onto convergence speed")
PyPlot.tight_layout()
