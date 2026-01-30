using FibonacciChain

N = 16
model = AnyonModel(IsingAnyon(), N; pbc=true)
H_sparse = anyon_ham_sparse(model)
energy, states = Arpack.eigs(H_sparse, nev=4, which=:SR)

GS = states[:,1]
σ = states[:,2]
ϵ = states[:,3]
∂σ = states[:,4]

corr_σσ = [abs(anyon_correlation(model, GS, σ, i, j)) for i in 1:N, j in 1:N]



formula(L) = -4/π*L - π/(6L)
Llis = collect(8:2:20)
E_data = formula.(Llis) 
inverseL = collect(0.00:0.0001: 0.015625) 
theoE = -4/π .- π/(6) .* inverseL
fig = plot(inverseL, theoE, xlabel=L"1/L^2", linewidth=2,
    label=false, ylabel=L"E_0^{+}/L", grid=false, 
    legend_background_color=nothing, legend_foreground_color=nothing, colorbar=false)
scatter!(fig, 1 ./(Llis.^2), E_data ./ Llis, label=false, markersize=4)