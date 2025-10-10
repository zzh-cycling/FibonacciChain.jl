using Arpack
using FibonacciChain

function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

function Boundarypost_selection_scaling(N)
    γlis = vcat(collect(0.0:0.05:0.95), [0.99, 0.999], 1.0)
    τlis = vcat(atanh.(vcat(collect(0.0:0.05:0.95), [0.99, 0.999])), 1e3)
    final_state_p_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    final_state_m_tau_lis = Vector{Vector{Float64}}(undef, length(τlis))
    final_sequence_p_tau_lis = Vector{Vector{Int64}}(undef, length(τlis))
    final_sequence_m_tau_lis = Vector{Vector{Int64}}(undef, length(τlis))
    total_weight_p_tau_lis = Vector{Float64}(undef, length(τlis))
    total_weight_m_tau_lis = Vector{Float64}(undef, length(τlis))
    EE_tau_lis_p=Vector{Vector{Float64}}(undef, length(τlis))
    EE_tau_lis_m=Vector{Vector{Float64}}(undef, length(τlis))
    FE_tau_lis_p = Vector{Float64}(undef, length(τlis))
    FE_tau_lis_m = Vector{Float64}(undef, length(τlis))

    @time energy, states = eigs(anyon_ham_sparse(N), nev=1, which=:SR)
    antiGS= states[:, 1]
    measurement_sites = collect(2:2:N)

    for (idx, τ) in enumerate(τlis)
        
        @time final_state_p, final_sequence_p, total_weight_p = boundary_post_selection(N, τ, antiGS, measurement_sites, 0)
        final_state_m, final_sequence_m, total_weight_m = boundary_post_selection(N, τ, antiGS, measurement_sites, 1)
        myprint(stdout, "N = $N, τ = $τ")

        EE_lis_p = eelis_Fibo_state(N, final_state_p)
        EE_lis_m = eelis_Fibo_state(N, final_state_m)
        FE_p = -log(total_weight_p)
        FE_m = -log(total_weight_m)
        
        EE_tau_lis_p[idx] = EE_lis_p
        EE_tau_lis_m[idx] = EE_lis_m
        FE_tau_lis_p[idx] = FE_p
        FE_tau_lis_m[idx] = FE_m
        final_state_p_tau_lis[idx] = final_state_p
        final_state_m_tau_lis[idx] = final_state_m
        final_sequence_p_tau_lis[idx] = final_sequence_p
        final_sequence_m_tau_lis[idx] = final_sequence_m
        total_weight_p_tau_lis[idx] = total_weight_p
        total_weight_m_tau_lis[idx] = total_weight_m
    end

    save("./exm/data/Boundary_postselection_N$(N).jld", "final_state_p_tau_lis", final_state_p_tau_lis, "final_state_m_tau_lis", final_state_m_tau_lis, "final_sequence_p_tau_lis", final_sequence_p_tau_lis, "final_sequence_m_tau_lis", final_sequence_m_tau_lis, "total_weight_p_tau_lis", total_weight_p_tau_lis, "total_weight_m_tau_lis", total_weight_m_tau_lis)
    save("./exm/data/Boundary_postselection__EE_FE_N$(N).jld", "EE_tau_lis_p", EE_tau_lis_p, "EE_tau_lis_m", EE_tau_lis_m, "FE_tau_lis_p", FE_tau_lis_p, "FE_tau_lis_m", FE_tau_lis_m)
    
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