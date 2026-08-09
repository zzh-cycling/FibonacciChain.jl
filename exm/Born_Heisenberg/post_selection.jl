using Distributed
using FibonacciChain
using JLD2
using Plots
using LinearAlgebra
include("../FitEntEntScal.jl")

if nprocs() == 1
    addprocs()
end

const BORN_HEISENBERG_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BORN_HEISENBERG_CONFIG)

@everywhere begin
    using FibonacciChain
    using Plots
    using Arpack
    using LinearAlgebra
    include("../FitEntEntScal.jl")

    function get_dynamics_params(ind, Δ)
        if ind == 1
            cfg = Dict(2.0 => (200, 2, 30))
            t, step, start = get(cfg, Δ, (100, 2, 30))
        elseif ind == 7
            cfg = Dict(2.0 => (25, 2, 6))
            t, step, start = get(cfg, Δ, (10, 2, 6))
        end
        inds = collect(1:step:2t)
        avg_range = start:2
        return t, inds, avg_range
    end

    function run_task_exact((Δ, N), ind)
        try
            τ = τlis[ind]
            model = heisenberg_model(N, Δ)
            st = ones(length(anyon_basis(model)))
            st ./= norm(st)
            t = get_dynamics_params(ind, Δ)[1]
            measure_config = MeasureConfig(
                τ = τ,
                t₂ = t*N,
                mode = :sample,
                cutoff = 1e-12,
                maxdim = 1000,
                verbose = false,
            )
            # All-false post-selection: rows = layers_per_period(model) * t₂ = 2tN layers,
            # cols = _samples_per_layer(model) = N÷2 (one sample per measured bond).
            samples = BitMatrix(zeros(Int8, 2t*N, N÷2))
            measure_outcome = bulk_evolution(model, st, measure_config, samples)
            F = measure_outcome.free_energys
            S_tlis = measure_outcome.entanglement_entropys
            final_st = measure_outcome.state
            Slis = anyon_eelis(model, final_st)
            mctable = Dict(8 => 2, 10 => 3)
            mc = get(mctable, N, 4)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut = mc, pbc = true)
            path = "exm/data/Heisenberg/Post_selection/eescaling_figs/gammaind$(ind)/L$(N)/Heisenberg_EntScal_Δ=$(round(Δ, digits=4))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (
                Δ = Δ,
                N = N,
                c = cent,
                c_err = cent_err,
                Slis = Slis,
                S_t = S_tlis,
                status = :success,
                error = nothing,
            )
        catch e
            return (Δ = Δ, N = N, status = :failed, error = e)
        end
    end
end

if length(ARGS) == 0
    println("Please provide L and τ index as command line arguments.")
else
    N = parse(Int64, ARGS[1])
    τ_idx = parse(Int64, ARGS[2])

    # Refined around the Δ = ±1 critical points of the XXZ chain.
    Δlis_sel =
        unique!(sort(vcat(collect(-2.0:0.1:2.0), collect(-1.02:0.02:-0.78), collect(0.78:0.02:1.02))))
    tasks = [(Δ, N) for Δ in Δlis_sel] |> vec

    results = pmap(tasks) do params
        run_task_exact(params, τ_idx)
    end

    failed_tasks = [(r.Δ, r.N, r.error) for r in results if r.status != :success]
    failed_count = length(failed_tasks)
    success_tasks =
        [(r.Δ, r.N, r.c, r.c_err, r.Slis, r.S_t) for r in results if r.status == :success]
    success_count = length(success_tasks)
    println("number of successful tasks: $success_count")
    println("number of failed tasks: $failed_count")
    println("failed tasks: $failed_tasks")


    cc_ensemble = collect([r[3] for r in success_tasks])
    cc_err_ensemble = collect([r[4] for r in success_tasks])
    Slis_ensemble = collect([r[5] for r in success_tasks])
    S_t_ensemble = collect([r[6] for r in success_tasks])
    save(
        "exm/data/Heisenberg/Post_selection/gammaind$(τ_idx)/ee_ensemble_L$(N).jld2",
        "Δlis",
        Δlis_sel,
        "cc_ensemble",
        cc_ensemble,
        "cc_err_ensemble",
        cc_err_ensemble,
        "Slis_ensemble",
        Slis_ensemble,
        "S_t_ensemble",
        S_t_ensemble,
    )
end
