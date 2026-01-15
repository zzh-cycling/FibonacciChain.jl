using Distributed
using FibonacciChain
using JLD2
using Plots
include("exm/FitEntEntScal.jl")

# 启动多进程（可根据实际情况调整进程数）
if nprocs() == 1
    addprocs()
end

@everywhere begin
    using FibonacciChain
    using Plots
    include("exm/FitEntEntScal.jl")
    Llis = collect(12:2:18)
    λlis = unique!(sort(vcat(collect(0.0:0.1:3), collect(0.4:0.02:0.6))))
    tasks = [(λ, N) for λ in λlis, N in Llis]
    tasks = vec(tasks)  
    function run_task_mps(task)
        try
            λ, N = task
            τ = atanh(1/√2)
            model = AnyonModel(OBFAnyon(), N, λ = λ, pbc=true)
            ψ, sites = initial_mps(N)
            t = 10N
            measure_config = MeasureConfig(τ=τ, t₂=t, mode=:sample, cutoff=1e-12, maxdim=1000, verbose=false)
            samples = BitMatrix(zeros(Int8, 8t, N))
            measure_outcome = bulk_evolution(model, sites, ψ, measure_config, samples)
            statelis, F = measure_outcome.states, measure_outcome.free_energys
            final_st = statelis[end]
            Slis = anyon_eelis(model, final_st)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/figs/OBF_EntScal_λ=$(round(λ, digits=3))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, c=NaN, c_err=NaN, status=:failed, error=e)
        end
    end

    function run_task_ed(task)
        try
            λ, N = task
            τ = atanh(1/√2)
            model = AnyonModel(OBFAnyon(), N, λ = λ, pbc=true)
            st = zeros(length(anyon_basis(model)))
            st[1] = 1.0
            t = 10N
            measure_config = MeasureConfig(τ=τ, t₂=t, mode=:sample, cutoff=1e-12, maxdim=1000, verbose=false)
            samples = BitMatrix(zeros(Int8, 8t, N))
            measure_outcome = bulk_evolution(model, st, measure_config, samples)
            statelis, F = measure_outcome.states, measure_outcome.free_energys
            final_st = statelis[end]
            Slis = anyon_eelis(model, final_st)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/figs/OBF_EntScal_λ=$(round(λ, digits=3))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, status=:success, error=nothing)
        catch e
            return (λ=λ, N=N, c=NaN, c_err=NaN, status=:failed, error=e)
        end
    end
end

results = pmap(run_task_ed, tasks)


clis       = [r.c      for r in results]
cerrorlis  = [r.c_err  for r in results]
cc_ensemble = reshape(clis, (length(λlis), length(Llis)))
cc_err_ensemble = reshape(cerrorlis, (length(λlis), length(Llis)))

failed_tasks = [(r.λ, r.N, r.c, r.c_err, r.error)
                for r in results if r.status != :success]

success_count = count(r -> r.status == :success, results)
failed_count  = length(failed_tasks)

println("number of successful tasks: $success_count")
println("number of failed tasks: $failed_count")
println("failed tasks: $failed_tasks")

save("exm/data/cc_ensemble.jld2", "λlis", λlis, "Llis", Llis, "cc_ensemble", cc_ensemble, "cc_err_ensemble", cc_err_ensemble)