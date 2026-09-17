using Distributed
using FibonacciChain
using JLD2
using LinearAlgebra

if nprocs() == 1
    addprocs()
end

const BORN_HEISENBERG_CONFIG = joinpath(@__DIR__, "config.jl")
@everywhere include($BORN_HEISENBERG_CONFIG)

@everywhere begin
    using FibonacciChain
    using Arpack
    using LinearAlgebra
    using JLD2

    function get_dynamics_params(ind, Δ)
        if ind == 1
            if Δ < -0.0
                cfg = Dict(-1.0 => (40, 2, 10))
                t, step, start = get(cfg, Δ, (40, 2, 10))
            else
                cfg = Dict(2.0 => (20, 2, 3))
                t, step, start = get(cfg, Δ, (15, 2, 3))
            end
            
        elseif ind == 2
            # γ = 0.1: weak measurement
            if Δ < -0.0
                cfg = Dict(-1.0 => (6, 2, 1))
                t, step, start = get(cfg, Δ, (6, 2, 1))
            else
                cfg = Dict(2.0 => (2, 2, 1)) 
                t, step, start = get(cfg, Δ, (2, 2, 1))
            end
        else
            cfg = Dict(2.0 => (2, 2, 1)) 
            t, step, start = get(cfg, Δ, (2, 2, 1))
        end
        inds = collect(1:step:2t)
        avg_range = start:2
        return t, inds, avg_range
    end
# When we choose AFM scenarios, the post-selected dynamics is imaginary-time evolution under the Floquet operator U(τ) = e^{∓τH_odd}e^{∓τH_even}; the steady state is the ground state of an effective Hamiltonian H_eff(τ) that interpolates between the XXZ chain (γ→0) and a sum of bond projectors (γ→1), which is a RVB state, the S(l) is ln2 or 2ln2. BCH gives H_eff = H_XXZ + (τ/2)[H_odd, H_even] + O(τ²), and the corrections are staggered under one-site translation — the brickwork pattern explicitly injects bond dimerization into H_eff. The dimerization operator has scaling dimension x_d = K with K = π/[2(π − arccos Δ)], so when we tune the γ, whether there is phase transition is Δ-dependent:
# • Δ > −1/√2 (K < 2): crossover, not a transition. Dimerization is relevant, so strictly any γ > 0 opens a gap m ~ τ^{1/(2−K)} — the LL exists only at γ = 0. But the gap is unobservably small at small τ, so on finite systems you see a smooth drift from c ≈ 1 to a dimerized profile.
# • Δ < −1/√2 (K > 2): a genuine transition at finite γ_c(Δ) < 1. Dimerization is irrelevant, so the LL is stable for a finite range of γ; but the γ=1 endpoint is provably an exact VBS, so a KT-type transition must occur in between. At Δ = −0.9 the system still looks critical at γ = 0.99. γ_c(Δ) is squeezed exponentially close to 1 as Δ → −1, because the per-bond singlet/aligned splitting 2J(1+Δ) vanishes there and projection requires τ ≳ 1/(1+Δ). 
# for Δ < −1 the AFM projective fixed point is the GHZ cat, S(l) = ln2 for every cut
    function run_task_exact((Δ, N), ind, sign::Bool)
        try
            τ = τlis[ind]
            model = heisenberg_model(N, Δ)
            basis = anyon_basis(model)
            st = ones(length(basis))
            if sign && τ >= 1e2
                # γ = 1 AFM branch: the uniform state has exactly zero overlap with
                # the singlet sector on any bond, so the first projective measurement
                # annihilates it. Seed the two Néel configurations to give the
                # post-selected trajectory nonzero weight (for Δ < -1 they are killed
                # instead and the uniform component carries the evolution).
                T = eltype(basis)
                neel_a = sum(1 << (2k) for k in 0:(N÷2-1))      # …0101
                neel_b = sum(1 << (2k + 1) for k in 0:(N÷2-1))  # …1010
                st[searchsortedfirst(basis, T(neel_a))] += 1.0
                st[searchsortedfirst(basis, T(neel_b))] += 1.0
            end
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
            # All-true post-selection: rows = layers_per_period(model) * t₂ = 2tN layers,
            # cols = _samples_per_layer(model) = N÷2 (one sample per measured bond).
            samples = sign ? BitMatrix(ones(Int8, 2t*N, N÷2)) : BitMatrix(zeros(Int8, 2t*N, N÷2)) 
            measure_outcome = bulk_evolution(model, st, measure_config, samples)
            F = measure_outcome.free_energys
            S_tlis = measure_outcome.entanglement_entropys
            final_st = measure_outcome.state
            Slis = anyon_eelis(model, final_st)
           
            branch = sign ? "AFM" : "FM"
            path = "exm/data/Heisenberg/Post_selection/$branch/gammaind$(ind)/L$(N)/Δ=$(round(Δ, digits=4))_N=$(N).jld2"
            mkpath(dirname(path))
            save(path,
                 "Slis", Slis, 
                 "S_tlis", S_tlis,
                 "F", F,
            )
            return (
                Slis = Slis,
                S_t = S_tlis,
                F = F,
                status = :success,
                error = nothing,
            )
        catch e
            return (Δ = Δ, N = N , status = :failed, error = e)
        end
    end
end

if length(ARGS) == 0
    println("Please provide L and τ index as command line arguments.")
else
    N = parse(Int64, ARGS[1])
    τ_idx = parse(Int64, ARGS[2])
    
    # Refined around the Δ = ±1 critical points of the XXZ chain.
    tasks = [(Δ, N) for Δ in Δlis] |> vec

    results = pmap(tasks) do params
        run_task_exact(params, τ_idx, true)
    end

    failed_tasks = [(r.Δ, r.N, r.error) for r in results if r.status != :success]
    failed_count = length(failed_tasks)
    success_tasks =
        [(r.Slis, r.S_t, r.F) for r in results if r.status == :success]
    success_count = length(success_tasks)
    println("number of successful tasks: $success_count")
    println("number of failed tasks: $failed_count")
    println("failed tasks: $failed_tasks")
end
