# Run from the repository root:
# julia --project=. exm/benchmark_mps_issue36.jl [repeats=3] [truncate_every_events=1]
# Run the same script with --project=/path/to/baseline for a before/after comparison.
using FibonacciChain, ITensorMPS, Random, LinearAlgebra, Statistics

function benchmark_mps_issue36(repeats=3, stride=1)
    repeats > 0 || throw(ArgumentError("repeats must be positive"))
    stride > 0 || throw(ArgumentError("truncate_every_events must be positive"))
    BLAS.set_num_threads(1)
    function run_case(L, chi, seed; periods=L == 8 ? 80 : 2L)
        model = AnyonModel(FibonacciAnyon(), L; pbc=true)
        ψ, sites = initial_mps(L)
        config = MeasureConfig(τ=L == 8 ? atanh(1/sqrt(2)) : atanh(0.95),
            mode=:Born, t₂=periods, rng=MersenneTwister(seed), cutoff=1e-12,
            maxdim=chi, truncate_every_events=stride)
        return bulk_evolution(model, sites, ψ, config)
    end
    # Compile the evolution path before measuring; setup and gate caching remain timed.
    run_case(8, 64, 1; periods=2)
    println("Julia $VERSION; Julia threads=$(Threads.nthreads()); BLAS threads=$(BLAS.get_num_threads())")
    for (L, chi) in ((8, 64), (16, 32), (16, 64), (32, 64))
        times, bytes = Float64[], Int[]
        for seed in 1:repeats
            GC.gc()
            result = @timed run_case(L, chi, seed)
            push!(times, result.time)
            push!(bytes, result.bytes)
            println((; L, chi, seed, stride, seconds=result.time,
                allocated_bytes=result.bytes, final_maxlink=maxlinkdim(result.value.state)))
            flush(stdout)
        end
        println((; L, chi, stride, median_seconds=median(times),
            median_allocated_bytes=median(bytes)))
        flush(stdout)
    end
end

benchmark_mps_issue36(isempty(ARGS) ? 3 : parse(Int, ARGS[1]),
    length(ARGS) < 2 ? 1 : parse(Int, ARGS[2]))
