include("../Meszaros_Problems/ConvertProblem.jl")

module Compare

using OSSDPTypes, JLD, Converter
export SolverResult, updateResults!,loadMeszarosData,getMeszarosDim,meszarosFilenames, printStatus, SolverMetrics, ProblemMetrics



  mutable struct SolverResult
    iter::Array{Int64,1}
    status::Array{Symbol,1}
    objVal::Array{Float64,1}
    x::Array{Array{Float64}}
    runTime::Array{Float64,1}
    numProblems::Int64
    problemDim::Array{Int64,2}
    problemName::Array{String,1}
    problemType::String
    solverName::String
    solverSettings::Array{Float64,1}
    timeStamp::String
    ind::Int64
    scalingON::Bool
    adaptionON::Bool
    objTrue::Float64
    xTrue::Array{Array{Float64}}

     #constructor
    function SolverResult(numProblems::Int64, problemType::String, solverName::String,timeStamp::String,solverSettings,scalingON,adaptionON)
    iter = zeros(Int64,numProblems)
    status = Array{Symbol}(numProblems)
    status[1:numProblems] = :empty
    objVal = zeros(numProblems)
    x = Array{Array{Float64}}(numProblems)
    runTime = zeros(numProblems)
    problemDim = zeros(Int64,numProblems,3)
    problemName = Array{String}(numProblems)
    problemName[1:numProblems] = "-"
    ind = 0
    objTrue = 0.
    xTrue = Array{Array{Float64}}(numProblems)

    if solverSettings == 0.
        settings = [solverSettings]
    else
        settings = [solverSettings.rho;solverSettings.sigma;solverSettings.alpha;solverSettings.scaling;solverSettings.eps_abs;solverSettings.eps_rel]
    end
    new(iter,status,objVal,x,runTime,numProblems,problemDim,problemName,problemType,solverName,settings,timeStamp,ind,scalingON,adaptionON,objTrue,xTrue)
    end
  end

   function Base.show(io::IO, obj::SolverResult)
    println(io,">>> Solver Results for: $(obj.solverName)\nTimeStamp: $(obj.timeStamp)\nProblem Type: $(obj.problemType)\nNumber of results:$(obj.ind)\nscaling ON: $(obj.scalingON)\nadaption ON: $(obj.adaptionON)" )
  end

  # define structs to hold metric data for each problem
  mutable struct SolverMetrics
    solverName::String
    adaptionON::Bool
    scalingON::Bool
    meanIterAll::Float64
    meanIterSolved::Float64
    numSolved::Int64
    percSolved::Float64
  end

  mutable struct ProblemMetrics
    problemType::String
    combinedData::Bool
    numProblems::Int64
    numSolvers::Int64
    solverResults::Array{SolverMetrics}
    function ProblemMetrics(problemType,combinedData,numProblems,numSolvers)
      new(problemType,combinedData,numProblems,numSolvers,Array{SolverMetrics}(numSolvers))
    end
  end


  function loadMeszarosData(data,solver::String)


    if contains(solver,"OSSDP")
      P, q, r, A, b, K = Converter.convertProblem(data)
      return P,q,r,A,b,K
    elseif contains(solver,"OSQP")
      P = data["P"]
      A = data["A"]
      q = data["q"]
      u = data["u"]
      l = data["l"]
      r = data["r"]
      return P, q[:],r,A,l[:], u[:]
    end
  end

  function getMeszarosDim(data)
    A = data["A"]
    return [size(A,1);size(A,2);nnz(A)]
  end

  function meszarosFilenames(path::String)
    fileNames = []
    for f in filter(x -> endswith(x, ".jld"), readdir(path))
        f = split(f,".")[1]
        push!(fileNames,String(f))
    end
    # sort filenames by number of nnz (stored in problemData[:,4])
    readmeInfo = JLD.load(path*"../objVals.jld")
    problemData = readmeInfo["problemData"]
    sortedInd = sort!(collect(1:1:length(fileNames)), by=i->problemData[i,4])
    fileNames = fileNames[sortedInd]

    # filter some problems by name
    excludeProbs = ["BOYD1";"BOYD2";"CONT-200";"CONT-201";"CONT-300";"UBH1";"QAFIRO";"QADLITTL"]
    filter!(x->(!in(x,excludeProbs)),fileNames)
    return fileNames
  end

  function printStatus(iii,numProblems,problem,resData)
    println("$(iii)/$(numProblems): $(problem) completed.")
    for jjj=1:length(resData)
      r = resData[jjj]
      r.scalingON ? si = " - scaled" : si = " - unscaled"
      println(" "^6*"$(r.solverName)$(si): Iterations: $(r.iter[r.ind]), Cost:$(r.objVal[r.ind]), Status:$(r.status[r.ind]), Runtime: $(r.runTime[r.ind])")
    end
  end

  function updateResults!(fn::String,resData,resArr,pDims::Array{Int64},pName::String,r::Union{Float64, Int64},SAVE_ALWAYS::Bool)
    numSolvers = length(resArr)
    n = pDims[2]

    for i=1:numSolvers
      resObj = resData[i]
      resObj.ind+=1
      resObj.problemDim[resObj.ind,:] = pDims
      resObj.problemName[resObj.ind] = pName
      if contains(resObj.solverName,"QOCS")
        resObj.iter[resObj.ind] = resArr[i].iter
        resObj.objVal[resObj.ind] = resArr[i].cost + r
        resObj.x[resObj.ind] = resArr[i].x[1:n]
        resObj.runTime[resObj.ind] = resArr[i].solverTime
        resObj.status[resObj.ind] = resArr[i].status
      elseif contains(resObj.solverName,"OSQP")
        resObj.iter[resObj.ind] = resArr[i].info.iter
        resObj.objVal[resObj.ind] = resArr[i].info.obj_val + r
        resObj.x[resObj.ind] = resArr[i].x
        resObj.status[resObj.ind] = resArr[i].info.status
      else
        nothing
      end
    end

    if SAVE_ALWAYS
      JLD.save(fn, "resData", resData)
    end


  end


end #module