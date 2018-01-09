# Test script to have an easy example solved by SCS that can be used to compare OSSDP solutions against
workspace()
include("../ossdp.jl")

using Base.Test
using OSSDP
using JLD
# using PyPlot


# # Problem DATA
A1 = [1.0 0 1; 0 3 7; 1 7 5]
A2 = [0.0 2 8; 2 6 0; 8 0 4]
C = [1.0 2 3; 2 9 0; 3 0 7]
b1 = 11.0
b2 = 19.0

P = zeros(3,3)
q = vec(C)
A = [vec(A1)';vec(A2)']
b = [b1;b2]

# σArr = [1e2,1e3]
# ρArr = [1,10]
# αArr = 0.2


σArr = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3,1e4,1e5,1e6,1e7]
ρArr = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3,1e4,1e5,1e6,1e7]
αArr = collect(0.1:0.1:1.9)
data = zeros(length(σArr)*length(ρArr)*length(αArr),6)

# main loop
counter = 1
for σ in σArr
  for ρ in ρArr
    for α in αArr
    settings = sdpSettings(rho=ρ,sigma=σ,alpha=α,max_iter=500,verbose=false)
    res,dbg = solveSDP(P,q,A,b,settings)
    normResult = norm(reshape(res.x,3,3)-[1.05593 0.369185 0.868298; 0.369185 0.129079 0.303585; 0.868299 0.303585 0.714011],Inf)
    normCost = norm(res.cost-13.902255839911007)
    normS = norm(res.x-res.s,Inf)

    data[counter,:] = [σ,ρ,α,normResult,normCost,normS]


    counter+=1
    end
  end
end


# analyze results

# parameter combination for lowest deviation from X
rowNum1 = findmin(data[:,4])
println("Lowest deviation from X: σ=$(data[rowNum1[2],1]), ρ=$(data[rowNum1[2],2]), α=$(data[rowNum1[2],3]), err=$(rowNum1[1])")
# parameter combination for lowest deviation from true cost
rowNum2 = findmin(data[:,5])
println("Lowest deviation from cost: σ=$(data[rowNum2[2],1]), ρ=$(data[rowNum2[2],2]), α=$(data[rowNum2[2],3]), err=$(rowNum2[1])")
# parameter combination for lowest deviation between X and S
rowNum3 = findmin(data[:,6])
println("Lowest deviation between X and S: σ=$(data[rowNum3[2],1]), ρ=$(data[rowNum3[2],2]), α=$(data[rowNum3[2],3]), err=$(rowNum3[1])")


# save to JLD file
filename = "parameterResults.jld"
JLD.save(filename, "data", data)
println("\nParameter Test Data successfully saved in $(filename).\n")

