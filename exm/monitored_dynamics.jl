function monitored_dynamics(L::Int64, τ::Float64, D::Int64, sign::Symbol=:m)
    L = 10
    D = 300L
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)
    
    samples_num = 100

    average_EE_tlis= zeros(D) 
    for i in 1:samples_num
        @show i
        sample_measured_states, samples, sample_weights = Bulkmeasure(L, τ, st, D) 
        EElis = zeros(D)
        for j in 1:D
            state_t = sample_measured_states[j]
            EE = eelis_Fibo_state(L, state_t)[5]
            EElis[j] = EE
        end
        final_state = sample_measured_states[end]
        average_EElis .+= eelis_Fibo_state(L, final_state)
        average_EE_tlis .+= EElis
    end

    average_EElis ./= samples_num
    average_EE_tlis ./= samples_num
    
    return average_EElis, average_EE_tlis
end

function test_post_selection(L::Int64, τ::Float64, D::Int64, sign::Symbol=:m)
    pbc = true
    st=zeros(length(Fibonacci_basis(L)))
    st[1] = 1.0
    average_EElis=zeros(L-1)

    EE_tlis = zeros(D)
    @time sample_measured_states, samples, sample_weights = Bulkpost_selection(L, τ, st, D, sign, pbc)
    for j in 1:D
        @show j
        state_t = sample_measured_states[j]
        EE = eelis_Fibo_state(L, state_t)[5]
        EE_tlis[j] = EE
    end
    final_state = sample_measured_states[end]
    average_EElis = eelis_Fibo_state(L, final_state)

    
    return average_EElis, EE_tlis
end

L = 10
τ = 0.1
D = 80L
average_EElis, EE_tlis = test_post_selection(L, τ, D, :m)
fig1=plot(range(0, D/2, D), EE_tlis, label=false, xlabel=L"t", ylabel=L"S_{vN}", title="Post-selection 0")
# plot!([1, 30*10], 0.8933161189003952*[1,1],  c=:Gray, label=false, linestyle=:dash, linewidth=2)
cent, fig = fitCCEntEntScal(average_EElis, mincut=2, pbc=true)
display(fig)
# save("exm/data/monitored_dynamics.jld", "cent", cent, "fig", fig)
# save("exm/data/monitored_dynamics.jld", "average_EElis", average_EElis, "EE_tlis", EE_tlis)
# savefig(fig, "/Users/cycling/Documents/projects/NoisyFibonacciChain/figs/Bulk_measure/monitored_dynamics_N28_D15_eescaling.pdf")
# savefig(fig1, "/Users/cycling/Documents/projects/NoisyFibonacciChain/figs/Bulk_measure/monitored_dynamics_N28_D15_S.pdf")