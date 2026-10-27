# Run from the repository root: julia --project exm/native_fibonacci_pbc.jl
using FibonacciChain, LinearAlgebra, Random

model = AnyonModel(FibonacciAnyon(), 4; pbc=true)
for flux in (:trivial, :tau)
    # Dense references use different complex conventions in the τ sector.
    H = native_pbc_dense(model; flux, construction=:mpo)
    H_braid = native_pbc_dense(model; flux, construction=:braid)
    @assert eigvals(Hermitian(H)) ≈ eigvals(Hermitian(H_braid))
    ψ, energy = anyon_mps_gst(model; backend=:native_pbc, flux,
        maxdim=32, sweep_times=12, cutoff=1e-12, outputlevel=0, seed=51)
    @assert isapprox(energy, eigmin(Hermitian(H)); atol=2e-9)
    println("flux = $flux, ED = ", eigvals(Hermitian(H)))
    println("DMRG = $energy")

    config = MeasureConfig(τ=0.8, t₂=2, mode=:Born,
        rng=MersenneTwister(36), maxdim=32, cutoff=1e-12)
    result = bulk_evolution(model, ψ.sites, ψ, config)
    @assert isapprox(norm(result.state), 1; atol=1e-10)
    println("measurement record = ", result.samples)
end
