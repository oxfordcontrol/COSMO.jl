
workspace()
using JLD

nameValDict = Dict()
# extract solution for objective value from README file
f = open("./MAT_FILES/README.txt")
lines = readlines(f)
close(f)
lines = lines[75:212]

data = Array{Any}(length(lines),5)
iii = 1
for ln in lines
  s = split(ln)
  prob = s[1]
  val = parse(Float64, s[7])
  nameValDict[prob] = val
  data[iii,:] = [s[1] parse(Int64, s[2]) parse(Int64, s[3]) parse(Int64, s[4]) parse(Float64, s[7])]
  iii += 1
end




  # save to JLD file
  JLD.save("./MAT_FILES/objVals.jld", "nameValDict",nameValDict,"problemData",data)
  println("Successfully saved obj name dictionary!")

