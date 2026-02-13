using Distributed
using Statistics
using FibonacciChain
using JLD
using Random
    
println("requested workers: ", nworkers())
println("total procs:       ", nprocs())
@everywhere println("host: ", gethostname(), "  pid: ", getpid())

@everywhere begin
    using FibonacciChain
    using JLD
    using Random
    using Statistics
    include("../FitEntEntScal.jl")

    binary_distribution(p, rng) = rand(rng) < p ? 1 : 0

    function get_system_params(τ, L)
        cfg = Dict(
            atanh(0.1)  => (2500L, 1000, 1500L),
            atanh(0.2)  => (500L,  100, 250L),
            atanh(0.3)  => (120L,  48, 100L),
            atanh(0.4)  => (100L,  40, 80L),
            atanh(0.5)  => (80L,   32, 40L),
            atanh(0.6)  => (45L,   20, 30L),
            log(1 + √2) => (35L,   14, 20L),
            atanh(0.8)  => (25L,   10, 10L),
            atanh(0.9)  => (8L,    4, 4L),
            atanh(0.95) => (8L,    4, 4L),
            atanh(0.999)=> (5L,    2, 2L),
        )
        D, step, start = get(cfg, τ, (5L, 2, 2L))
        return D, collect(1:step:D), start:D-5
    end

    function ps_prob_evolution(params)
        L, τ, seed = params
        D, _, _ = get_system_params(τ, L)
        model = AnyonModel(FibonacciAnyon(), L; pbc=true)
        problis = collect(0.1:0.1:0.9)
        ee_plis = Vector{Vector{Float64}}(undef, length(problis))
        initial_state = zeros(length(anyon_basis(model)))
        initial_state[1] = 1.0
        gate_num = div(D*L, 2)

        for (idx, prob) in enumerate(problis)
            rng = MersenneTwister(seed)
            sample = BitMatrix(reshape([binary_distribution(prob, rng) for _ in 1:gate_num], D, div(L, 2)))
            config = MeasureConfig(τ=τ, mode=:sample, t₂=div(D,2))
            mo = bulk_evolution(model, initial_state, config, sample)
            ee_plis[idx] = anyon_eelis(model, mo.state)
        end
        
        mkpath("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)")
        save("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_sample$(seed).jld", 
             "seed", seed, "ee_plis", ee_plis, "problis", problis)
             
        return ee_plis
    end
    # function ps_prob_evolution_Ising(params)
        # Try to reproduce the outcome of `Entanglement Transition in the Projective Transverse Field Ising Model`, but fail, due to it's not measurement at everywhere.
    function process_ps_prob_evolution(L, τ)
        # Load the data
        D = get_system_params(τ, L)[1]
        samplelis = collect(1:10000)
        problis = collect(0.1:0.1:0.9)
        centlis = []
        seedlis = zeros(Int64, length(samplelis))
        # Process the data
        ee_problis = zeros(L-1, length(problis))
        for (i, sample) in enumerate(samplelis)
            ee_plis, seed = load("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_sample$(sample).jld", "ee_plis", "seed")
            ee_problis += hcat(ee_plis...)
            seedlis[i] = seed
        end
    
        ee_problis ./= length(samplelis)
        
    
        for i in eachindex(problis)
            push!(centlis, fitCCEntEntScal(vec(ee_problis[:, i]), mincut=2, pbc=true)[1])
        end
    
        save("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_cent.jld", "centlis", centlis, "seedlis", seedlis)
        return centlis, seedlis
    end
    
    function process_data()
        Llis = collect(8:2:20)
        problis = collect(0.1:0.1:0.9)
        ixs = [1, 3, 4, 7, 9, 10, 12]
        for (inds, τ) in enumerate(τlis[ixs])
            cent_Lplis = zeros(length(Llis), length(problis))
            cent_stderrlis = zeros(length(Llis), length(problis))
            for (id, L) in enumerate(Llis)
                D, _, _ = get_system_params(τ, L)
                @show (L, τ)
                centlis= load("exm/data/Bulk_measure/ps_prob_evolution/L$(L)/τ$(τ)/L$(L)_D$(div(D,L))_τ$(τ)_cent.jld", "centlis")
                cent_Lplis[id, :] = [i[1] for i in centlis]
                cent_stderrlis[id, :] = [i[2] for i in centlis]
            end
            save("exm/data/Bulk_measure/ps_prob_evolution/centlis_L$(Llis[1])$(Llis[end])_τ$(τ).jld", "cent_Lplis", cent_Lplis, "cent_stderrlis", cent_stderrlis)
        end
    end
end
    
  

γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.707, 0.8, 0.9, 0.95, 0.999, 1]
τlis = atanh.(γlis)
τlis[end] = 1000.0  # Last value is for γ=1
τlis[findfirst(γlis .== 0.707)] = log(1 + √2) 
seed_interval_lis = collect(1:100:2000)

if length(ARGS) == 0
    println("No arguments provided.")
    println("Usage: julia -p N ps_prob.jl L_start L_end τ_idx seed_start seed_end")
else
    L_start = parse(Int64, ARGS[1])
    L_end = parse(Int64, ARGS[2])
    inds = parse(Int64, ARGS[3])
    seed_start = parse(Int64, ARGS[4])
    seed_end = parse(Int64, ARGS[5])
    
    τ = τlis[inds]
    Llis = collect(L_start:2:L_end)
    seeds = collect(seed_start:seed_end)
    
    # 生成所有 (L, τ, seed) 组合
    tasks = [(L, τ, seed) for L in Llis for seed in seeds]
    println("Running $(length(tasks)) tasks: L=$Llis, τ=$τ, seeds=$seed_start:$seed_end on $(nprocs()) workers")
    
    results = @time pmap(tasks) do params
        try
            ps_prob_evolution(params)
            return (params, :success, nothing)
        catch e
            return (params, :failed, e)
        end
    end
    
    # statistics on results
    succeeded = filter(r -> r[2] == :success, results)
    failed = filter(r -> r[2] == :failed, results)
    println("\n=== Statistics ===")
    println("Success: $(length(succeeded)) / $(length(tasks))")
    println("Failed: $(length(failed)) / $(length(tasks))")
    
    if !isempty(failed)
        failed_params = [r[1] for r in failed]
        println("Failed params: $failed_params")
    end
end