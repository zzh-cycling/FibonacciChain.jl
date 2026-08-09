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

    function run_task_GS_ed((Δ, N))
        try
            model = heisenberg_model(N, Δ)
            H = anyon_ham_sparse(model)
            energy, states = Arpack.eigs(H, nev = 1, which = :SR)
            GS = states[:, 1]
            Slis = anyon_eelis(model, GS)
            (cent, cent_err), fig = fitCCEntEntScal(Slis, mincut = 4, pbc = true)
            path = "exm/data/Heisenberg/GS/eescaling_figs/Heisenberg_Δ=$(round(Δ, digits=3))_N=$(N).pdf"
            mkpath(dirname(path))
            savefig(fig, path)
            return (
                Δ = Δ,
                N = N,
                c = cent,
                c_err = cent_err,
                Slis = Slis,
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

    # Refined around the Δ = ±1 critical points of the XXZ chain; J = +1 (AFM side).
    Δlis_sel =
        unique!(sort(vcat(collect(-2.0:0.1:2.0), collect(-1.02:0.02:-0.78), collect(0.78:0.02:1.02))))
    tasks = [(Δ, N) for Δ in Δlis_sel] |> vec

    results = pmap(tasks) do params
        run_task_GS_ed(params)
    end

    failed_tasks = [(r.Δ, r.N, r.error) for r in results if r.status != :success]
    failed_count = length(failed_tasks)
    success_tasks =
        [(r.Δ, r.N, r.c, r.c_err, r.Slis) for r in results if r.status == :success]
    success_count = length(success_tasks)
    println("number of successful tasks: $success_count")
    println("number of failed tasks: $failed_count")
    println("failed tasks: $failed_tasks")


    cc_ensemble = collect([r[3] for r in success_tasks])
    cc_err_ensemble = collect([r[4] for r in success_tasks])
    Slis_ensemble = collect([r[5] for r in success_tasks])
    mkpath("exm/data/Heisenberg/GS/L$(N)")
    save(
        "exm/data/Heisenberg/GS/L$(N)/GS_cc_ensemble.jld2",
        "Δlis",
        Δlis_sel,
        "cc_ensemble",
        cc_ensemble,
        "cc_err_ensemble",
        cc_err_ensemble,
        "Slis_ensemble",
        Slis_ensemble,
    )

end
