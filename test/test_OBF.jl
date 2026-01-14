using FibonacciChain
using Test
using LinearAlgebra
using BitBasis
using Arpack
using Random 
using LsqFit

@testset "measure_matrix" begin
    
    # measuring XZZ (X_i Z_{i+1} Z_{i+2})
    τ = 1.0
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    model_XZZ = AnyonModel(OBFAnyon(), N, pbc=true, measure_operator=:XZZ)
    # At idx=1: X_1 Z_2 Z_3
    expected_matrix = cstτ * I(8) + coef * σx ⊗ σz ⊗ σz
    Mpobc = FibonacciChain.measure_matrix(model_XZZ, τ, 1, false)
    @test Mpobc ≈ expected_matrix
    
    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ * I(8) + coef * σx ⊗ σz ⊗ σz
    Mmobc = FibonacciChain.measure_matrix(model_XZZ, τ, 1, true)
    @test Mmobc ≈ expected_matrix
    @test Mpobc^2 + Mmobc^2 ≈ I(8)

    # measuring ZZX (Z_i Z_{i+1} X_{i+2})
    coef = sinh(τ/2) / √(2cosh(τ))
    model_ZZX = AnyonModel(OBFAnyon(), N, pbc=true, measure_operator=:ZZX)
    # At idx=1: Z_1 Z_2 X_3
    expected_matrix = cstτ * I(8) + coef * σz ⊗ σz ⊗ σx
    Mpobc = FibonacciChain.measure_matrix(model_ZZX, τ, 1, false)
    @test Mpobc ≈ expected_matrix
    
    coef = -sinh(τ/2) / √(2cosh(τ))
    expected_matrix = cstτ * I(8) + coef * σz ⊗ σz ⊗ σx
    Mmobc = FibonacciChain.measure_matrix(model_ZZX, τ, 1, true)
    @test Mmobc ≈ expected_matrix
    @test Mpobc^2 + Mmobc^2 ≈ I(8)

    # Test XZZ/ZZX with larger system and different index
    N4 = 4
    model_XZZ4 = AnyonModel(OBFAnyon(), N4, pbc=true, measure_operator=:XZZ)
    model_ZZX4 = AnyonModel(OBFAnyon(), N4, pbc=true, measure_operator=:ZZX)
    # At idx=2: X_2 Z_3 Z_4
    expected_XZZ = cstτ * I(16) + sinh(τ/2) / √(2cosh(τ)) * I(2) ⊗ σx ⊗ σz ⊗ σz
    @test FibonacciChain.measure_matrix(model_XZZ4, τ, 2, false) ≈ expected_XZZ
    # At idx=2: Z_2 Z_3 X_4
    expected_ZZX = cstτ * I(16) + sinh(τ/2) / √(2cosh(τ)) * I(2) ⊗ σz ⊗ σz ⊗ σx
    @test FibonacciChain.measure_matrix(model_ZZX4, τ, 2, false) ≈ expected_ZZX
end

@testset "measuremap OBFAnyon" begin
    # Test OBFAnyon measuremap with :XZZ, :ZZX, and :OBF operators
    N = 4
    τ = 1.0
    sign = false
    
    # Common coefficients
    cstτ = cosh(τ/2) / √(2cosh(τ))
    coef = sinh(τ/2) / √(2cosh(τ))
    
    # Test :XZZ operator
    model_xzz = AnyonModel(OBFAnyon(), N; pbc=false, measure_operator=:XZZ)
    state = ones(2^N)
    idx = 1
    output_xzz = measuremap(model_xzz, τ, state, idx, sign)
    
    # XZZ flips site 1, with coefficient depending on ZZ eigenvalue of sites 2,3
    # For basis states where sites 2,3 have same spin: zz_eigen = +1
    # For basis states where sites 2,3 have different spin: zz_eigen = -1
    @test length(output_xzz) == 2^N
    @test sum(output_xzz) > 0  # State should be non-trivial
    
    # Test :ZZX operator
    model_zzx = AnyonModel(OBFAnyon(), N; pbc=false, measure_operator=:ZZX)
    output_zzx = measuremap(model_zzx, τ, state, idx, sign)
    
    # ZZX flips site 3 (idx+2), with coefficient depending on ZZ eigenvalue of sites 1,2
    @test length(output_zzx) == 2^N
    @test sum(output_zzx) > 0
    
    # Test :OBF operator (should be equivalent to XZZ followed by ZZX)
    model_obf = AnyonModel(OBFAnyon(), N; pbc=false, measure_operator=:OBF)
    output_obf = measuremap(model_obf, τ, state, idx, sign)
    
    # OBF = XZZ then ZZX, so output should match sequential application
    state_after_xzz = measuremap(model_xzz, τ, state, idx, sign)
    state_after_both = measuremap(model_zzx, τ, state_after_xzz, idx, sign)
    @test output_obf ≈ state_after_both
    
    # Test with different initial state
    state2 = collect(1.0:2^N)
    output_obf2 = measuremap(model_obf, τ, state2, idx, sign)
    state2_after_xzz = measuremap(model_xzz, τ, state2, idx, sign)
    state2_after_both = measuremap(model_zzx, τ, state2_after_xzz, idx, sign)
    @test output_obf2 ≈ state2_after_both
    
    # Test with sign = true
    sign = true
    output_obf_sign = measuremap(model_obf, τ, state, idx, sign)
    state_xzz_sign = measuremap(model_xzz, τ, state, idx, sign)
    state_both_sign = measuremap(model_zzx, τ, state_xzz_sign, idx, sign)
    @test output_obf_sign ≈ state_both_sign
    
    # Test PBC
    model_obf_pbc = AnyonModel(OBFAnyon(), N; pbc=true, measure_operator=:OBF)
    model_xzz_pbc = AnyonModel(OBFAnyon(), N; pbc=true, measure_operator=:XZZ)
    model_zzx_pbc = AnyonModel(OBFAnyon(), N; pbc=true, measure_operator=:ZZX)
    
    sign = false
    output_obf_pbc = measuremap(model_obf_pbc, τ, state, idx, sign)
    state_xzz_pbc = measuremap(model_xzz_pbc, τ, state, idx, sign)
    state_both_pbc = measuremap(model_zzx_pbc, τ, state_xzz_pbc, idx, sign)
    @test output_obf_pbc ≈ state_both_pbc
    
    # Test τ → ∞ limit (projective measurement)
    τ_large = 1e3
    model_obf_proj = AnyonModel(OBFAnyon(), N; pbc=false, measure_operator=:OBF)
    model_xzz_proj = AnyonModel(OBFAnyon(), N; pbc=false, measure_operator=:XZZ)
    model_zzx_proj = AnyonModel(OBFAnyon(), N; pbc=false, measure_operator=:ZZX)
    
    output_obf_proj = measuremap(model_obf_proj, τ_large, state, idx, sign)
    state_xzz_proj = measuremap(model_xzz_proj, τ_large, state, idx, sign)
    state_both_proj = measuremap(model_zzx_proj, τ_large, state_xzz_proj, idx, sign)
    @test output_obf_proj ≈ state_both_proj
    
    # Verify OBF preserves normalization structure
    normalized_state = ones(2^N) / √(2^N)
    output_norm = measuremap(model_obf, τ, normalized_state, idx, sign)
    @test norm(output_norm) > 0  # Should produce non-zero output
end

function fitCCEntEntScal(
    SvN_list::Vector{Float64};
    err::Vector{Float64}=0.0SvN_list,
    mincut::Int=1,
    pbc::Bool=false)

    # log of chord length / 6 for open boundary
    logChord(l, L) = @. log(sin(π * l /L))/6
    
    L = length(SvN_list) + 1

    # fit scaling
    lm(x,p) = @. p[1] * x + p[2]
    xdata = logChord([1:L-1;],L); #log.(sin.(π .* [1:L-1;] ./L))./6
    fit = curve_fit(lm, xdata[mincut:L-mincut], SvN_list[mincut:L-mincut], [0.5, 0.0])
    fitparam = fit.param
    cent = fitparam[1]
    cent_err = stderror(fit)[1]
    if pbc
        cent /= 2.0
        cent_err/= 2.0
    end
    println("cent ± cent_err is $(cent) ± $(cent_err)")

    return (cent, cent_err)
end

@testset "boundary_evolution" begin
    
end

@testset "bulk_evolution" begin
    
end

@testset "central_charge" begin
    N = 12
    τ = atanh(1/√2) # critical point for IsingX
    model = AnyonModel(OBFAnyon(), N, λ = 0.856, pbc=true)

    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0
    t = 6N
    samples = BitMatrix(zeros(Int8, 8t, N))
    measure_config = MeasureConfig(τ=τ, t₂=t, mode=:sample)

    measure_outcome = bulk_evolution(model, st, measure_config, samples)
    generated_statelis, F = measure_outcome.states, measure_outcome.free_energys
    final_st = generated_statelis[end]
    EE = anyon_eelis(model, final_st)
    @test fitCCEntEntScal(EE, mincut=2, pbc=true)[1][1] ≈ 0.8 atol=1e-1

end 