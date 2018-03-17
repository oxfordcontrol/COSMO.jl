# Test routine to compare performance for a number of SDP problems (max Eigenvalue problems)
workspace()
include("../../src/Solver.jl")
include("../solverComparison/Compare.jl")

using OSQP, OSSDP, Base.Test, Compare

rng = MersenneTwister(123734)

nn = 15
timestamp = Dates.format(now(), "yyddmm_HH-MM")
dataFile = "./SC_" * timestamp * "MaxEigSDP.jld"
problemType = "SDP-MaxEig"
setOFF = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=false)
setON = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=false)
setOFFAd = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=true)
setONAd = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=1500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-5,eps_rel=1e-5,adaptive_rho=true)
sr1 = SolverResult(nn, problemType,"QOCS",timestamp,setOFF,false)
sr2 = SolverResult(nn, problemType,"QOCS",timestamp,setON,true)
sr3 = SolverResult(nn, problemType,"QOCS-Ad",timestamp,setOFF,false)
sr4 = SolverResult(nn, problemType,"QOCS-Ad",timestamp,setON,true)
resData = [sr1;sr2;sr3;sr4]
ws = 0
for iii =1:1:nn
 # generate symmetric test matrix A
  nA = rand(rng,2:35)
  A = Symmetric(randn(rng,nA,nA))

  # define max eigenvalue problem as SDP

  # create augmented matrices
  Aa = -vec(eye(nA))
  ba = -vec(A)
  q = [1]
  P = spzeros(1,1)
  # define cone membership
  K = Cone(0,0,[],[nA^2])
  pDims = [size(Aa,1);size(Aa,2);size(Aa,1)*size(Aa,2)]


  # Solve with OSSDP
  resOSSDP_unscaled,nothing = OSSDP.solve(P,q,Aa,ba,K,setOFF);
  print("\n.")
  resOSSDP_scaled,nothing = OSSDP.solve(P,q,Aa,ba,K,setON);
  print(".")
   resOSSDP_unscaledAd,nothing = OSSDP.solve(P,q,Aa,ba,K,setOFFAd);
  print("\n.")
  resOSSDP_scaledAd,nothing = OSSDP.solve(P,q,Aa,ba,K,setONAd);
  print(".")

  println("Correct max eigenvalue: $(maximum(eig(A)[1]))")

  resArray = [resOSSDP_unscaled;resOSSDP_scaled;resOSSDP_unscaledAd;resOSSDP_scaledAd]
  updateResults!(dataFile,resData,resArray,pDims,"SDP-$(iii)",0.,true)
  printStatus(iii,nn,"SDP-$(iii)",resData)
end


