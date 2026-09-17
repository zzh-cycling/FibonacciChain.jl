# Run from the repository root:
# julia --project=. exm/benchmark_mps_issue36.jl [repeats=3] [truncate_every_events=1] [Fibonacci|Ising|OBF]
# Run the same script with --project=/path/to/baseline for a before/after comparison.
using FibonacciChain, ITensorMPS, Random, LinearAlgebra, Statistics

function benchmark_mps_issue36(repeats=3, stride=1, model_type=:Fibonacci)
    repeats > 0 || throw(ArgumentError("repeats must be positive"))
    stride > 0 || throw(ArgumentError("truncate_every_events must be positive"))
    model_type in (:Fibonacci, :Ising, :OBF) || throw(ArgumentError("unknown model: $model_type"))
    BLAS.set_num_threads(1)
    fibonacci = model_type == :Fibonacci
    function run_case(L, chi, seed; periods=fibonacci ? (L == 8 ? 80 : 2L) : 12)
        model = fibonacci ? AnyonModel(FibonacciAnyon(), L; pbc=true) :
            AnyonModel(SpinHalf(), L; model_type, pbc=true, λ=0.3)
        ψ, sites = initial_mps(L)
        τ = fibonacci ? (L == 8 ? atanh(1/sqrt(2)) : atanh(0.95)) : 0.8
        config = MeasureConfig(τ=τ,
            mode=:Born, t₂=periods, rng=MersenneTwister(seed), cutoff=1e-12,
            maxdim=chi, truncate_every_events=stride)
        return bulk_evolution(model, sites, ψ, config)
    end
    # Compile the evolution path before measuring; setup and gate caching remain timed.
    run_case(fibonacci ? 8 : 6, 64, 1; periods=2)
    println("Julia $VERSION; Julia threads=$(Threads.nthreads()); BLAS threads=$(BLAS.get_num_threads())")
    println("Model: $model_type")
    cases = fibonacci ? ((8, 64), (16, 32), (16, 64), (32, 64)) : ((6, 16), (24, 32), (24, 64))
    for (L, chi) in cases
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
    length(ARGS) < 2 ? 1 : parse(Int, ARGS[2]),
    length(ARGS) < 3 ? :Fibonacci : Symbol(ARGS[3]))
