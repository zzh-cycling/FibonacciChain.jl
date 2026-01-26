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
      function run_task_GS_ed((λ, N))
        try
            if λ >= 10.0
                model = AnyonModel(OBFAnyon(), N; λI=0.0, pbc=true)
            else
                model = AnyonModel(OBFAnyon(), N; λ=λ, pbc=true)
            end
            H = anyon_ham_sparse(model)
            energy, states = Arpack.eigs(H, nev=1, which=:SR)
            GS = states[:, 1]
            Slis = anyon_eelis(model, GS)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut=4, pbc=true)
            path = "exm/data/OBF/GS/eescaling_figs/OBF_λ=$(round(λ, digits=3))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (λ=λ, N=N, c=cent, c_err=cent_err, Slis=Slis, status=:success, error=nothing)
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
    
    λlis = unique!(sort(vcat(collect(0.0:0.1:2), collect(0.816:0.04:1.02),[11.0])))
    tasks  = [(λ, N) for λ in λlis] |> vec
    
    results = pmap(tasks) do params
        run_task_GS_ed(params)
    end

    failed_tasks = [(r.λ, r.N, r.error)
                    for r in results if r.status != :success]
    failed_count  = length(failed_tasks)
    success_tasks = [(r.λ, r.N, r.c, r.c_err, r.Slis)
                    for r in results if r.status == :success]
    success_count = length(success_tasks)
    println("number of successful tasks: $success_count")
    println("number of failed tasks: $failed_count")
    println("failed tasks: $failed_tasks")
    

    cc_ensemble = collect([r[3]  for r in success_tasks])
    cc_err_ensemble  = collect([r[4] for r in success_tasks])
    Slis_ensemble = collect([r[5] for r in success_tasks])
    mkpath("exm/data/OBF/GS/L$(N)")
    save("exm/data/OBF/GS/L$(N)/GS_cc_ensemble.jld2", "λlis", λlis, "cc_ensemble", cc_ensemble, "cc_err_ensemble", cc_err_ensemble, "Slis_ensemble", Slis_ensemble)

end