using Distributed
using FibonacciChain
using JLD
using Statistics

# Start parallel workers if not already running
if nprocs() == 1
    addprocs()
end

@everywhere begin
    using FibonacciChain
    using JLD
    using Statistics
    
    function post_selection(L::Int64, τ::Float64, D::Int64, sign::Bool=true)
        pbc = true
        model = AnyonModel(FibonacciAnyon(), L; pbc=pbc)
        st=zeros(length(anyon_basis(model)))
        st[1] = 1.0
        average_EElis=zeros(L-1)

        sample = sign ? BitMatrix(ones(Bool, D, length(2:2:L))) : BitMatrix(zeros(Bool, D, length(2:2:L)))
        config = MeasureConfig(τ=τ, mode=:sample, t₂=div(D,2))
        mo = bulk_evolution(model, st, config, sample)
        sample_measured_states, sample, sample_free_energy = mo.states, mo.samples, mo.free_energys
        EE_tlis = [ee(anyon_rdm(model, collect(1:div(L,2)), state_t)) for state_t in sample_measured_states]
        final_state = sample_measured_states[end]
        average_EElis = anyon_eelis(model, final_state)

        
        return average_EElis, EE_tlis, sample_free_energy
    end

    function get_system_params(τ, L)
        table = Dict(
                atanh(0.1)  => (300L, 1000, 1500L),
                atanh(0.2)  => (60L,  100, 250L),
                atanh(0.3)  => (25L,  48, 100L),
                atanh(0.4)  => (20L,  40, 80L),
                atanh(0.5)  => (20L,   32, 40L),
                atanh(0.6)  => (15L,   20, 30L),
                log(1 + √2) => (10L,   14, 20L),
                atanh(0.8)  => (8L,   10, 10L),
                atanh(0.9)  => (5L,    4, 4L),
                atanh(0.95) => (5L,    4, 4L),
            )
        D, step, start = get(table, τ, (3L, 2, 2L))
        inds = collect(1:step:D)
        avg_range = start:D-5
        return D, inds, avg_range
    end
end

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1, and atanh(1/√2) = log(1 + √2)

function run_task(task)
    i, L = task
    τ = τlis[i]
    D = get_system_params(τ, L)[1]
    sign = true
    
    println("Processing: i=$i, L=$L, τ=$τ, D=$D")
    
    average_EElis, EE_tlis, sample_free_energy = post_selection(L, τ, D, sign)
    
    # Save results
    filepath = "exm/data/post_selection$(sign)/τ$(τ)/L$(L)_D$(div(D,L)).jld"
    mkpath(dirname(filepath))
    save(filepath, "average_EElis", average_EElis, "EE_tlis", EE_tlis, "sample_free_energy", sample_free_energy)
    
    println("Completed: i=$i, L=$L")
    return (i=i, L=L, τ=τ, status="completed")
end

if length(ARGS) == 0
    println("No arguments provided.")
else
    # Build task list: all combinations of (i, L)
    tasks = [(i, L) for i in 1:12 for L in [22, 24]]
    
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