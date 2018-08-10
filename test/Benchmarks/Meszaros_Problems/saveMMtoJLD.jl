# Script used to convert problem data in .mat files into .jld files

using JLD, MAT
dirPath = "./bart_meszaros_data/"

problems = []
for f in filter(x -> endswith(x, ".mat"), readdir(dirPath))
    f = split(f,".")[1]
    push!(problems,String(f))
end
readmeInfo = JLD.load("./MAT_FILES/objVals.jld")
problemData = readmeInfo["problemData"]
sortedInd = sort!(collect(1:1:length(problems)), by=i->problemData[i,4])
problems = problems[sortedInd]

numProblems = length(problems)



nameValDict = JLD.load("./MAT_FILES/objVals.jld")["nameValDict"]



iii = 1
for iii=1:length(problems)
  problem = problems[iii]
  data = MAT.matread("$(dirPath)"*"$(problem).mat")

    P = data["P"]
    A = data["A"]
    newA = sparse(full(A))
    q = data["q"]
    u = data["u"]
    l = data["l"]
    r = data["r"]
    objVal = nameValDict[lowercase(problem)]

    JLD.save("./JLD_Files/"*problem*".jld", "P",P,"A",newA,"q",q,"r",r,"u",u,"l",l,"objVal",objVal)
    println("$(iii)/$(length(problems))")
end
