using Distributed
using FibonacciChain
using JLD
using Statistics

# Start parallel workers if not already running
if nprocs() == 1
    addprocs()
end

const BULK_MEASURE_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BULK_MEASURE_CONFIG)

@everywhere begin
    using FibonacciChain
    using JLD
    using Statistics

    function post_selection(L::Int64, τ::Float64, D::Int64, sign::Bool = true)
        pbc = true
        model = AnyonModel(FibonacciAnyon(), L; pbc = pbc)
        st=zeros(length(anyon_basis(model)))
        st[1] = 1.0
        average_EElis=zeros(L-1)

        sample =
            sign ? BitMatrix(ones(Bool, D, length(2:2:L))) :
            BitMatrix(zeros(Bool, D, length(2:2:L)))
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = div(D, 2))
        mo = bulk_evolution(model, st, config, sample)
        sample, sample_free_energy = mo.samples, mo.free_energys
        EE_tlis = mo.entanglement_entropys
        final_state = mo.state
        average_EElis = anyon_eelis(model, final_state)


        return average_EElis, EE_tlis, sample_free_energy
    end

    function post_selection_mps(L::Int64, τ::Float64, χ::Int64, sign::Bool = true)
        model = fib_model(L)
        t, _, _ = get_mps_ps(τ, L, sign)
        average_EElis=zeros(L-1)
        ψ, sites = initial_mps(L)
            
        sample =
            sign ? BitMatrix(ones(Bool, 2t, length(2:2:L))) :
            BitMatrix(zeros(Bool, 2t, length(2:2:L)))
        config = MeasureConfig(τ = τ, mode = :sample, t₂ = t, cutoff = 1e-12, maxdim = χ,)

        try
            @time mps_mo = bulk_evolution(model, sites, ψ, config, sample)
            sample, sample_free_energy = mps_mo.samples, mps_mo.free_energys
            EE_tlis = mps_mo.entanglement_entropys
    
            # Compute final state EE profile
            average_EElis = anyon_eelis(model, mps_mo.state)
            
            branch = sign ? "AFM" : "FM"
            filepath = "exm/data/post_selection$(branch)/τ$(τ)/L$(L)_t$(div(t, L))_χ$(χ).jld"
            mkpath(dirname(filepath))
            save(
                filepath,
                "average_EElis",
                average_EElis,
                "EE_tlis",
                EE_tlis,
                "sample_free_energy",
                sample_free_energy,
            )
            return (L, τ, χ, :success, nothing)
        catch e
            return (L, τ, χ, :failed, e)
        end
    end

    function get_mps_ps(τ, L, sign)
        if sign == true
            table = Dict(
                atanh(0.1) => (50L, 1000, 7L),
                atanh(0.2) => (40L, 100, 12L),
                atanh(0.3) => (30L, 48, 5L),
                atanh(0.4) => (20L, 40, 4L),
                atanh(0.5) => (15L, 32, 3L),
                atanh(0.6) => (12L, 20, 2L),
                log(1 + √2) => (8L, 14, L),
                atanh(0.8) => (6L, 10, L),
                atanh(0.9) => (2L, 4, L),
                atanh(0.95) => (L, 4, L),
            )
            # No idea why it will drop after L+δt, maybe because exact dynamics preserve the Fibonacci constraint algebraically, but ordinary SVD truncation does not explicitly enforce it
            t, step, start = get(table, τ, (L, 2, L))
            inds = collect(1:step:div(t, 2))
            avg_range = start:(t-4)
            return t, inds, avg_range
        else
            table = Dict(
                atanh(0.1) => (100L, 1000, 15L),
                atanh(0.2) => (40L, 100, 5L),
                atanh(0.3) => (30L, 48, 4L),
                atanh(0.4) => (25L, 40, 3L),
                atanh(0.5) => (20L, 32, 3L),
                atanh(0.6) => (12L, 20, 2L),
                log(1 + √2) => (8L, 14, 2L),
                atanh(0.8) => (4L, 10, 1L),
                atanh(0.9) => (2L, 4, L),
                atanh(0.95) => (L, 4, L),
            )
            t, step, start = get(table, τ, (L, 2, 2L))
            inds = collect(1:step:t)
            avg_range = start:(t-4)
            return t, inds, avg_range
        end
    end

    function get_system_params_ps(τ, L, sign)
        if sign == true
            table = Dict(
                atanh(0.1) => (2500L, 1000, 750L),
                atanh(0.2) => (500L, 100, 120L),
                atanh(0.3) => (120L, 48, 50L),
                atanh(0.4) => (100L, 40, 40L),
                atanh(0.5) => (80L, 32, 20L),
                atanh(0.6) => (45L, 20, 15L),
                log(1 + √2) => (35L, 14, 10L),
                atanh(0.8) => (25L, 10, 5L),
                atanh(0.9) => (8L, 4, 2L),
                atanh(0.95) => (8L, 4, 2L),
                atanh(0.999) => (5L, 2, 1L),
            )
            D, step, start = get(table, τ, (5L, 2, L))
            inds = collect(1:step:div(D, 2))
            avg_range = start:(div(D, 2)-5)
            return D, inds, avg_range
        else
            table = Dict(
                atanh(0.1) => (300L, 1000, 1500L),
                atanh(0.2) => (60L, 100, 250L),
                atanh(0.3) => (25L, 48, 100L),
                atanh(0.4) => (20L, 40, 80L),
                atanh(0.5) => (20L, 32, 40L),
                atanh(0.6) => (15L, 20, 30L),
                log(1 + √2) => (10L, 14, 20L),
                atanh(0.8) => (8L, 10, 10L),
                atanh(0.9) => (5L, 4, 4L),
                atanh(0.95) => (5L, 4, 4L),
            )
            D, step, start = get(table, τ, (3L, 2, 2L))
            inds = collect(1:step:D)
            avg_range = start:(D-5)
            return D, inds, avg_range
        end
    end

    function run_task(task)
        i, L = task
        τ = τlis[i]
        sign = true
        D = get_system_params_ps(τ, L, sign)[1]

        println("Processing: i=$i, L=$L, τ=$τ, D=$D")

        average_EElis, EE_tlis, sample_free_energy = post_selection(L, τ, D, sign)

        # Save results
        branch = sign ? "AFM" : "FM"
        filepath = "exm/data/post_selection$(branch)/τ$(τ)/L$(L)_D$(div(D,L)).jld"
        mkpath(dirname(filepath))
        save(
            filepath,
            "average_EElis",
            average_EElis,
            "EE_tlis",
            EE_tlis,
            "sample_free_energy",
            sample_free_energy,
        )

        println("Completed: i=$i, L=$L")
        return (i = i, L = L, τ = τ, status = "completed")
    end
end


if length(ARGS) == 0
    println("No arguments provided.")
else
    # Build task list: all combinations of (i, L)
    tasks = [(i, L) for i = 1:12 for L in [22, 24]]

    println("Total tasks: $(length(tasks))")
    println("Number of workers: $(nworkers())")

    # Run tasks in parallel
    results = pmap(run_task, tasks)

    println("\nAll tasks completed!")
    println("Results summary:")
    for r in results
        println("  i=$(r.i), L=$(r.L), τ=$(r.τ) - $(r.status)")
    end
end
