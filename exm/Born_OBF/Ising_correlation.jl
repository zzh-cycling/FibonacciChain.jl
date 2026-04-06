using FibonacciChain
using Arpack
using LinearAlgebra
using Plots
using LaTeXStrings

function vacuum_energy(Llis)
    E0lis = Float64[]
    for L in Llis
        model = AnyonModel(IsingAnyon(), L; pbc=true)
        H_sparse = anyon_ham_sparse(model)
        energy, states = Arpack.eigs(H_sparse, nev=1, which=:SR)
        push!(E0lis, real(energy[1]))
    end
    return E0lis
end

function vacuum_energy_casimir(Llis)
    E0lis = vacuum_energy(Llis)
    
    formula(L) = -4/π*L - π/(6L)
    E_data = formula.(Llis) 
    inverseL = collect(0.00:0.0001: 0.015625) 
    theoE = -4/π .- π/(6) .* inverseL
    fig = plot(1 ./Llis, theoE, xlabel=L"1/L^2", linewidth=2,
        label=false, ylabel=L"E_0^{+}/L", grid=false, 
        legend_background_color=nothing, legend_foreground_color=nothing, colorbar=false)
    scatter!(fig, 1 ./(Llis.^2), E_data ./ Llis, label=false, markersize=4)

    return fig
end

function spectrum_momentum(N)
    model = AnyonModel(IsingAnyon(), N; pbc=true)
    klis = collect(0:N-1)
    energy_spectrum_momentum = zeros(N, 10)

    for k in klis
        H_sparse_k = anyon_ham_sparse(model, k)
        energy, states = Arpack.eigs(H_sparse_k, nev=10, which=:SR)
        energy_spectrum_momentum[k+1, :] = real(energy)
    end
    
    return klis, energy_spectrum_momentum
end    

function conformal_tower(N)
    klis, energy_spectrum_momentum = spectrum_momentum(N)
    E0 = minimum(energy_spectrum_momentum)
    v = 2.0  # velocity for Ising CFT

    fig = plot(title="Conformal Tower for Ising CFT", xlabel="Momentum k", ylabel="(E - E0) * (N / 2πv)", legend=false)
    for n in 1:size(energy_spectrum_momentum, 2)
        scaled_energies = (energy_spectrum_momentum[:, n] .- E0) .* (N / (2π * v))
        scatter!(fig, klis, scaled_energies, markersize=4)
    end

    return fig
end

function ising_correlation(N)
    model = AnyonModel(IsingAnyon(), N; pbc=true)
    H_sparse = anyon_ham_sparse(model)
    energy, states = Arpack.eigs(H_sparse, nev=1, which=:SR)
    GS = states[:, 1]
    corr_σσ = [abs(anyon_correlation(model, GS, σ, i, j)) for i in 1:N, j in 1:N]

    return corr_σσ
end

function correlation_scaling_dimension(N)
    corr_σσ = ising_correlation(N)
    rlis = collect(1:div(N,2))
    C_r = [corr_σσ[1, r+1] for r in rlis]

    fig = plot(rlis, C_r, seriestype=:scatter, xlabel=L"r", ylabel=L"C(r) = \langle \sigma_0 \sigma_r \rangle", title="Spin-Spin Correlation Function Scaling", label=false)
    fit_func(r, A, Δ) = A * r .^ (-2 * Δ)
    p_fit = curve_fit(fit_func, rlis, C_r, [1.0, 0.125])
    A_fit, Δ_fit = p_fit.param

    r_fit = collect(1:0.1:div(N,2))
    C_fit = fit_func(r_fit, A_fit, Δ_fit)
    plot!(fig, r_fit, C_fit, label="Fit: Δ=$(round(Δ_fit, digits=4))", linewidth=2, color=:red)

    return (Δ=Δ_fit, fig=fig)
end