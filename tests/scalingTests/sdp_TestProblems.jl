# Test routine to compare scaling for a number of max eigenvalue problems (SDPs)

workspace()
include("../../src/Solver.jl")

using Base.Test, JLD
using OSSDP, OSSDPTypes

nn = 25
rng = MersenneTwister(7232)
resIter = zeros(nn,6)
resCost = zeros(nn,6)
resDim = zeros(nn)

for iii = 1:nn
    # generate symmetric test matrix A
    r = rand(rng,2:50)
    A = Symmetric(randn(rng,r,r))

    # solve the dual problem
    ca = -vec(A)
    Aa = vec(eye(r))'
    ba = [1.]
    K = Cone(0,0,[],[r^2])
    Pa = zeros(r^2,r^2)


    # Solve with OSSDP
    settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
    resOSSDP_unscaled,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
    settings = OSSDPSettings(rho=1.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
    resOSSDP_scaled,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
    settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
    resOSSDP_unscaled100,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
    settings = OSSDPSettings(rho=1000.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 0,eps_abs = 1e-3,eps_rel=1e-4)
    resOSSDP_unscaled1000,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
    settings = OSSDPSettings(rho=100.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
    resOSSDP_scaled100,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)
    settings = OSSDPSettings(rho=1000.0,sigma=1.0,alpha=1.6,max_iter=2500,verbose=false,checkTermination=1,scaling = 10,eps_abs = 1e-3,eps_rel=1e-4)
    resOSSDP_scaled1000,nothing = OSSDP.solve(Pa,ca,Aa,ba,K,settings)


    resCost[iii,:] = [resOSSDP_unscaled.cost resOSSDP_scaled.cost resOSSDP_unscaled100.cost resOSSDP_scaled100.cost resOSSDP_unscaled1000.cost resOSSDP_scaled1000.cost]
    resIter[iii,:] = [resOSSDP_unscaled.iter resOSSDP_scaled.iter resOSSDP_unscaled100.iter resOSSDP_scaled100.iter resOSSDP_unscaled1000.iter resOSSDP_scaled1000.iter]
    resDim[iii] = r
    println("$(iii)/$(nn) completed!")
end


timestamp = Dates.format(now(), "yyddmm_HH-MM")
filename = timestamp * "_sdpScalingTest.jld"
JLD.save(filename, "resCost", resCost, "resIter",resIter,"resDim",resDim)
println(">>> Test Data successfully saved in $(filename).\n")
