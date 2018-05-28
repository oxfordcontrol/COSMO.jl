using PyPlot, JLD

fn = "qocs_vs_osqp_meszaros.jld"
data = JLD.load(fn)
n = length(data["iters"])

println("""QOCS success ratio: $(data["success"]/n)""")
println("""OSQP success ratio: $(data["successOSQP"]/n)""")

using PyPlot
figure(1)
ratio = data["iters"][2:end]./data["itersOSQP"][2:end]
semilogy(2:n, ratio, linestyle="None", marker=".")
idx = []
max_iter = []
for i in 1:length(ratio)
    if data["iters"][i+1] == 5000 || data["itersOSQP"][i+1] == 5000
        append!(idx, i + 1)
        append!(max_iter, ratio[i])
    end
end
semilogy(idx, max_iter, c="r", linestyle="None", marker=".")

title("Number of Iterations in the Maros-Meszaros Dataset")
# legend(["QOCS", "OSQP"])
xlabel("Problem index (sorted by nnz)")
ylabel("Iterations QOCS / Iterations OSQP")
savefig("iter.pdf")

figure(2)
semilogy(2:n, 1000*[data["times"][2:end] data["timesOSQP"][2:end]], linestyle="None", marker=".")
title("Average per iteration time for the Maros-Mescaros Dataset")
legend(["QOCS", "OSQP"])
xlabel("Problem index (sorted by nnz)")
ylabel("Time (ms)")
savefig("timings.png")
