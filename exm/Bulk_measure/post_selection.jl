using Distributed
using FibonacciChain
using JLD2
using Plots
using LinearAlgebra
include("../FitEntEntScal.jl")

if nprocs() == 1
    addprocs()
end

@everywhere begin
    using FibonacciChain
    using Plots
    using Arpack
    using LinearAlgebra
    include("../FitEntEntScal.jl")
    γlis = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/√2, 0.8, 0.9, 0.95, 0.999, 1]
    τlis = atanh.(γlis)
    τlis[end] = 1000.0
    τlis[findfirst(γlis .== 1/√2)] = log(1 + √2)
    function get_dynamics_params(ind, λ)
        if ind == 1
            cfg = Dict(
                0.936  => (150, 1000, 750),
                0.976 => (600,   14, 10),
                1.0   => (2000,   14, 10),
                1.016 => (1500,   14, 10),
                1.1  => (800,   14, 10),
                1.2 => (500,   14, 10),
                1.3 => (400,   14, 10),
                1.7 => (4000,   14, 10),
                11.0 => (4000,   14, 10),
            )
            t, step, start = get(cfg, λ, (100, 14, 60))
        elseif ind == 7
            if 0.7< λ <0.776
                t, step, start = (100, 14, 30)
            elseif λ == 0.776
                t, step, start = (150, 14, 30)
            elseif 0.776< λ < 0.9
                t, step, start = (200, 14, 60)
            elseif 0.9 <=λ < 0.95
                t, step, start = (50, 14, 60)
            elseif λ == 1.2
                t, step, start = (25, 14, 60)
            elseif λ == 11.0
                t, step, start = (25, 14, 60)
            else
                t, step, start = (10, 14, 6)
            end
        end
        inds = collect(1:step:14t)
        avg_range = start:14-5
        return t, inds, avg_range
    end
    function run_task_mps((λ, N), ind)
        try
            τ = τlis[ind]
            if λ >= 10.0
                model = AnyonModel(OBFAnyon(), N; λI=0.0, pbc=true)
            else
                model = AnyonModel(OBFAnyon(), N; λ=λ, pbc=true)
            end
            ψ, sites = initial_mps(N)
            t = get_dynamics_params(ind, λ)[1]
            measure_config = MeasureConfig(τ=τ, t₂=t*N, mode=:sample, cutoff=1e-12, maxdim=1000, verbose=false)
            samples = BitMatrix(zeros(Int8, 14t*N, N))
            measure_outcome = bulk_evolution(model, sites, ψ, measure_config, samples)
            statelis, F = measure_outcome.states, measure_outcome.free_energys
            S_tlis = [ee_mps(st, div(N, 2)) for st in statelis]
            final_st = statelis[end]
            Slis = anyon_eelis(model, final_st)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/data/OBF/Dynamics/eescaling_figs/gammaind$(ind)/OBF_EntScal_λ=$(round(λ, digits=4))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, Slis=Slis, S_t=S_tlis, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, status=:failed, error=e)
        end
    end
  
    function run_task_exact((λ, N), ind)
        try
            τ = τlis[ind]
            if λ >= 10.0
                model = AnyonModel(OBFAnyon(), N; λI=0.0, pbc=true)
            else
                model = AnyonModel(OBFAnyon(), N; λ=λ, pbc=true)
            end
            #  Using even parity state will converge faster.
            st = ones(length(anyon_basis(model)))
            st ./= norm(st)
            t = get_dynamics_params(ind, λ)[1]
            measure_config = MeasureConfig(τ=τ, t₂=t*N, mode=:sample, cutoff=1e-12, maxdim=1000, verbose=false)
            samples = BitMatrix(zeros(Int8, 14t*N, N))
            measure_outcome = bulk_evolution(model, st, measure_config, samples)
            statelis, F = measure_outcome.states, measure_outcome.free_energys
            S_tlis = [ee(anyon_rdm(model, collect(1:div(N, 2)), st)) for st in statelis]
            final_st = statelis[end]
            Slis = anyon_eelis(model, final_st)
            mctable = Dict(
                8 => 2,
                10 => 3,
            )
            mc = get(mctable, N, 4)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=mc, pbc=true)
            path = "exm/data/OBF/Dynamics/eescaling_figs/gammaind$(ind)/OBF_EntScal_λ=$(round(λ, digits=4))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, Slis=Slis, S_t=S_tlis, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, status=:failed, error=e)
        end
    end
end

if length(ARGS) == 0
    println("Please provide L and τ index as command line arguments.")
else
    N = parse(Int64, ARGS[1])
    τ_idx = parse(Int64, ARGS[2])
    
    if τ_idx== 1
        λlis = unique!(sort(vcat(collect(0.0:0.1:2), collect(0.816:0.04:1.02),[11.0])))
    elseif τ_idx== 7
        λlis = unique!(sort(vcat(collect(0.0:0.1:1.5), collect(0.616:0.02:1.02),[11.0])))
    end
    tasks  = [(λ, N) for λ in λlis] |> vec
    
    results = pmap(tasks) do params
        run_task_exact(params, τ_idx)
    end

    failed_tasks = [(r.λ, r.N, r.error)
                    for r in results if r.status != :success]
    failed_count  = length(failed_tasks)
    success_tasks = [(r.λ, r.N, r.c, r.c_err, r.Slis, r.S_t)
                    for r in results if r.status == :success]
    success_count = length(success_tasks)
    println("number of successful tasks: $success_count")
    println("number of failed tasks: $failed_count")
    println("failed tasks: $failed_tasks")
    

    cc_ensemble = collect([r[3]  for r in success_tasks])
    cc_err_ensemble  = collect([r[4] for r in success_tasks])
    Slis_ensemble = collect([r[5] for r in success_tasks])
    S_t_ensemble = collect([r[6] for r in success_tasks])
    mkpath("exm/data/OBF/Dynamics/gammaind$(τ_idx)/L$(N)")
    save("exm/data/OBF/Dynamics/gammaind$(τ_idx)/L$(N)/GS_cc_ensemble2.jld2", "λlis", λlis, "cc_ensemble", cc_ensemble, "cc_err_ensemble", cc_err_ensemble, "Slis_ensemble", Slis_ensemble, "S_t_ensemble", S_t_ensemble)

end
