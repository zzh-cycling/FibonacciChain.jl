# julia --project=. exm/benchmark_mps_backends.jl [L=8] [chi=32] [periods=4]
# Same fixed record and encoded entropy; timings exclude compilation.
using FibonacciChain, Random, LinearAlgebra
import MPSKit, TensorKit, ITensorMPS

function benchmark_mps_backends(L=8, chi=32, periods=4)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and >= 4"))
    chi > 0 && periods > 0 || throw(ArgumentError("chi and periods must be positive"))
    BLAS.set_num_threads(1)
    model = AnyonModel(FibonacciAnyon(), L)
    kit, ks = initial_mps(L; backend=:mpskit)
    old, os = initial_mps(L; backend=:itensor)
    born = MeasureConfig(τ=0.8, t₂=periods, mode=:Born, rng=MersenneTwister(36), maxdim=chi)
    record = bulk_evolution(model, os, old, born).samples
    config = MeasureConfig(τ=0.8, t₂=periods, mode=:sample, maxdim=chi)
    run(sites, state) = bulk_evolution(model, sites, state, config, record)
    run(ks, kit)
    run(os, old)
    results = Dict()
    for (name, sites, state) in ((:itensor, os, old), (:mpskit, ks, kit))
        GC.gc()
        result = @timed run(sites, state)
        results[name] = result.value
        ψ = result.value.state
        bond = name == :itensor ? ITensorMPS.maxlinkdim(ψ) :
            maximum(TensorKit.dim(MPSKit.right_virtualspace(ψ, i)) for i in 1:L-1)
        println((; backend=name, L, chi, periods, seconds=result.time,
            allocated_bytes=result.bytes, max_bond=bond, final_entropy=last(result.value.entanglement_entropys)))
    end
    overlap = abs(dot(results[:mpskit].state, mpskit_state(results[:itensor].state)))^2
    println((; fidelity=overlap, max_free_energy_difference=maximum(abs.(
        results[:mpskit].free_energys-results[:itensor].free_energys))))
    return nothing
end

benchmark_mps_backends(length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8,
    length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 32,
    length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 4)
