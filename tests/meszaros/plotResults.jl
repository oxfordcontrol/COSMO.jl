using PyPlot, JLD

fn = "qocs_vs_osqp_meszaros.jld"
data = JLD.load(fn)
n = length(data["filenames"])

println("""QOCS success ratio: $(data["success"]/n)""")
println("""OSQP success ratio: $(data["successOSQP"]n)""")

using PyPlot
figure(1)
semilogy(2:n, [data["iters"] data["itersOSQP"]], linestyle="None", marker=".")
title("Number of Iterations in the Maros-Mescaros Dataset")
legend(["QOCS", "OSQP"])
xlabel("Problem index (sorted by nnz)")
ylabel("Iterations")
savefig("iter.pdf")

figure(2)
semilogy(2:n, 1000*[data["times"] data["timesOSQP"]], linestyle="None", marker=".")
title("Average per iteration time for the Maros-Mescaros Dataset")
legend(["QOCS", "OSQP"])
xlabel("Problem index (sorted by nnz)")
ylabel("Time (ms)")
savefig("timings.pdf")
