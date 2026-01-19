using Distributed
using FibonacciChain
using JLD2
using Plots
include("../FitEntEntScal.jl")

if nprocs() == 1
    addprocs()
end

@everywhere begin
    using FibonacciChain
    using Plots
    using Arpack
    include("../FitEntEntScal.jl")
    γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
    τlis = atanh.(γlis)
    τlis[end] = 1000.0
    τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)
    function get_dynamics_params(ind, L)
        cfg = Dict(
            1  => (100L, 1000, 750L),
            7 => (7L,   14, 10L),
        )
        t, step, start = get(cfg, ind, (5L, 2, L))
        inds = collect(1:step:14t)
        avg_range = start:14-5
        return t, inds, avg_range
    end
    function run_task_mps((λ, N), ind)
        try
            τ = τlis[ind]
            model = AnyonModel(OBFAnyon(), N, λ = λ, pbc=true)
            ψ, sites = initial_mps(N)
            t = get_dynamics_params(ind, N)[1]
            measure_config = MeasureConfig(τ=τ, t₂=t, mode=:sample, cutoff=1e-12, maxdim=1000, verbose=false)
            samples = BitMatrix(zeros(Int8, 14t, N))
            measure_outcome = bulk_evolution(model, sites, ψ, measure_config, samples)
            statelis, F = measure_outcome.states, measure_outcome.free_energys
            final_st = statelis[end]
            Slis = anyon_eelis(model, final_st)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/data/OBF/Dynamics/eescaling_figs/gammaind$(ind)/OBF_EntScal_λ=$(round(λ, digits=4))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, c=NaN, c_err=NaN, status=:failed, error=e)
        end
    end
    function run_task_GS_ed((λ, N))
        try
            model = AnyonModel(OBFAnyon(), N, λ = λ, pbc=true)
            H = anyon_ham_sparse(model)
            energy, states = Arpack.eigs(H, nev=1, which=:SR)
            GS = states[:, 1]
            Slis = anyon_eelis(model, GS)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/data/OBF/GS/eescaling_figs/OBF_λ=$(round(λ, digits=3))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, c=NaN, c_err=NaN, status=:failed, error=e)
        end
    end
    function run_task_exact((λ, N), ind)
        try
            τ = τlis[ind]
            model = AnyonModel(OBFAnyon(), N, λ = λ, pbc=true)
            st = zeros(length(anyon_basis(model)))
            st[1] = 1.0
            t = get_dynamics_params(ind, N)[1]
            measure_config = MeasureConfig(τ=τ, t₂=t, mode=:sample, cutoff=1e-12, maxdim=1000, verbose=false)
            samples = BitMatrix(zeros(Int8, 14t, N))
            measure_outcome = bulk_evolution(model, st, measure_config, samples)
            statelis, F = measure_outcome.states, measure_outcome.free_energys
            final_st = statelis[end]
            Slis = anyon_eelis(model, final_st)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/data/OBF/Dynamics/eescaling_figs/gammaind$(ind)/OBF_EntScal_λ=$(round(λ, digits=4))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, c=NaN, c_err=NaN, status=:failed, error=e)
        end
    end
end

if length(ARGS) == 0
    println("Please provide L and τ index as command line arguments.")
else
    N = parse(Int64, ARGS[1])
    τ_idx = parse(Int64, ARGS[2])
    
    λlis = unique!(sort(vcat(collect(0.0:0.1:2), collect(0.816:0.04:1.02),[5.0])))
    tasks  = [(λ, N) for λ in λlis] |> vec
    
    results = pmap(tasks) do params
        run_task_exact(params, τ_idx)
    end

    cc_ensemble = collect([r.c  for r in results])
    cc_err_ensemble  = collect([r.c_err for r in results])

    failed_tasks = [(r.λ, r.N, r.c, r.c_err, r.error)
                    for r in results if r.status != :success]

    success_count = count(r -> r.status == :success, results)
    failed_count  = length(failed_tasks)

    println("number of successful tasks: $success_count")
    println("number of failed tasks: $failed_count")
    println("failed tasks: $failed_tasks")
    mkpath("exm/data/OBF/Dynamics/gammaind$(τ_idx)/L$(N)")
    save("exm/data/OBF/Dynamics/gammaind$(τ_idx)/L$(N)/GS_cc_ensemble.jld2", "λlis", λlis, "cc_ensemble", cc_ensemble, "cc_err_ensemble", cc_err_ensemble)
end