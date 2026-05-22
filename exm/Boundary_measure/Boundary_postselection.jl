using Arpack
using FibonacciChain

γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)

function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

function Boundarypost_selection_scaling(N)
    final_state_p_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    final_state_m_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    final_sequence_p_tau_lis = Vector{BitVector}(undef, length(τlis))
    final_sequence_m_tau_lis = Vector{BitVector}(undef, length(τlis))
    EE_tau_lis_p=Vector{Vector{Float64}}(undef, length(τlis))
    EE_tau_lis_m=Vector{Vector{Float64}}(undef, length(τlis))
    FE_tau_lis_p = Vector{Float64}(undef, length(τlis))
    FE_tau_lis_m = Vector{Float64}(undef, length(τlis))

    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    @time energy, states = eigs(anyon_ham_sparse(model), nev = 1, which = :SR)
    antiGS = states[:, 1]
    measurement_sites = collect(2:2:N)

    for (idx, τ) in enumerate(τlis)
        config_p = MeasureConfig(τ = τ, t₂ = 1, mode = :sample, enable_τ_eff = false)
        config_m = MeasureConfig(τ = τ, t₂ = 1, mode = :sample, enable_τ_eff = false)
        sample_p = BitVector(zeros(div(N, 2)))
        sample_m = BitVector(ones(div(N, 2)))
        # Note:  boundary_evolution with specific samples
        @time mo_p = boundary_evolution(model, antiGS, config_p, sample_p)
        state_p, final_sequence_p, free_energy_p = mo_p.state, mo_p.sample, mo_p.free_energy
        @time mo_m = boundary_evolution(model, antiGS, config_m, sample_m)
        state_m, final_sequence_m, free_energy_m = mo_m.state, mo_m.sample, mo_m.free_energy
        myprint(stdout, "N = $N, τ = $τ")

        EE_lis_p = anyon_eelis(model, state_p)
        EE_lis_m = anyon_eelis(model, state_m)

        EE_tau_lis_p[idx] = EE_lis_p
        EE_tau_lis_m[idx] = EE_lis_m
        FE_tau_lis_p[idx] = free_energy_p
        FE_tau_lis_m[idx] = free_energy_m
        final_state_p_tau_lis[idx] = state_p
        final_state_m_tau_lis[idx] = state_m
        final_sequence_p_tau_lis[idx] = final_sequence_p
        final_sequence_m_tau_lis[idx] = final_sequence_m
    end

    save(
        "./exm/data/Boundary_postselection_N$(N).jld",
        "final_state_p_tau_lis",
        final_state_p_tau_lis,
        "final_state_m_tau_lis",
        final_state_m_tau_lis,
        "final_sequence_p_tau_lis",
        final_sequence_p_tau_lis,
        "final_sequence_m_tau_lis",
        final_sequence_m_tau_lis,
    )
    save(
        "./exm/data/Boundary_postselection_EE_FE_N$(N).jld",
        "EE_tau_lis_p",
        EE_tau_lis_p,
        "EE_tau_lis_m",
        EE_tau_lis_m,
        "FE_tau_lis_p",
        FE_tau_lis_p,
        "FE_tau_lis_m",
        FE_tau_lis_m,
    )

end

if length(ARGS) == 0
    println("No arguments provided.")
else
    for arg in ARGS
        println("Received argument: $arg")
        N=parse(Int64, arg)
        Boundarypost_selection_scaling(N)
    end
end
