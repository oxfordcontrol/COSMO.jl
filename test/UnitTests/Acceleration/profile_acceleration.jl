using LinearAlgebra, Profile, COSMO, Random
using StatProfilerHTML, BenchmarkTools
Profile.init(n = 10^7, delay = 0.001)
total_iter = 2000
dim = 100000
mem = 10

w = zeros(dim)
w_prev = zeros(dim)

function run_profiling!(accelerator, total_iter, dim, mem, w, w_prev)
    rng = Random.MersenneTwister(1)

    for i = 1:total_iter
        # make new random inputs but reuse memory
        for i in eachindex(w)
            w[i] = 10. * rand(rng) .- 5.
            w_prev[i] = 10. * rand(rng) .- 5.
        end

        COSMO.update_history!(accelerator, w, w_prev, i)
        COSMO.accelerate!(w, w_prev, accelerator, i)
    end

end

# Accelerator with QR decomposition for Least Squares problem
accelerator = AndersonAccelerator{Float64, NoRegularizer, Type2{QRDecomp}, RestartedMemory}(dim, mem = mem)
accelerator.activated = true
@profilehtml run_profiling!(accelerator, total_iter, dim, mem, w, w_prev)

# Accelerator with Normal Equations for Least Squares problem
# Uncomment here to profile the accelerator with normal equations!
# accelerator2 = AndersonAccelerator{Float64, NoRegularizer, Type2{NormalEquations}, RestartedMemory}(dim, mem = mem)
# accelerator2.activated = true
# @profilehtml run_profiling!(accelerator2, total_iter, dim, mem, w, w_prev)




# Time both methods
# Uncomment here to time both methods against each other
# accelerator = AndersonAccelerator{Float64, NoRegularizer, Type2{QRDecomp}, RestartedMemory}(dim, mem = mem)
# accelerator.activated = true
# @btime run_profiling!(accelerator, total_iter, dim, mem, $w, $w_prev)

# accelerator2 = AndersonAccelerator{Float64, NoRegularizer, Type2{NormalEquations}, RestartedMemory}(dim, mem = mem)
# accelerator2.activated = true
# @btime run_profiling!(accelerator2, total_iter, dim, mem, $w, $w_prev)
