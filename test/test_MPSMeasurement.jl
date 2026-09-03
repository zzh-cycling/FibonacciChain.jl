using Test
using FibonacciChain
using ITensorMPS, ITensors
using LinearAlgebra
using Random
using Arpack
using BitBasis

@testset "initial_mps" begin
    N = 6
    # Test initial MPS generation
    ψ, sites = initial_mps(N)

    @test ψ isa MPS
    @test length(ψ) == N
    @test all(linkdims(ψ) .== 1)  # Initial state should have bond dimension 1

    # Test normalization, set to be 1
    @test abs(inner(ψ, ψ) - 1.0) < 1e-10
end

@testset "Fibonacci DMRG ground state and TCI central charge" begin
    L = 16
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    ψ, energy = anyon_mps_gst(
        model;
        sweep_times = 20,
        maxdim = 128,
        cutoff = 1e-12,
        outputlevel = 0,
        seed = 1234,
    )

    exact_energy = real(Arpack.eigs(anyon_ham_sparse(model); nev = 1, which = :SR)[1][1])
    violation = fibonacci_constraint_violation(ψ; pbc = true)
    @test ψ isa MPS
    @test length(ψ) == L
    @test inner(ψ, ψ) ≈ 1 atol = 1e-10
    @test energy ≈ exact_energy atol = 5e-4
    @test violation < 1e-6

    # For a periodic critical ground state,
    # S(l) = (c/3) log[(L/pi) sin(pi*l/L)] + constant.
    cuts = collect(4:(L - 4))
    chord = log.(L / π .* sin.(π .* cuts ./ L)) ./ 3
    entropy = anyon_eelis(model, ψ)[cuts]
    central_charge = (hcat(ones(length(cuts)), chord) \ entropy)[2]
    @test central_charge ≈ 0.7 atol = 0.04
end

@testset "Fibonacci dense/MPO Hamiltonian convention" begin
    L = 4
    sites = siteinds("Qubit", L)
    for pbc in (false, true), measure_operator in (:Antiferro, :Ferro)
        model = AnyonModel(
            FibonacciAnyon(),
            L;
            pbc = pbc,
            measure_operator = measure_operator,
        )
        basis = anyon_basis(model)
        product_states = [
            productMPS(
                sites,
                [string((state.buf >> (L - site)) & 1) for site = 1:L],
            ) for state in basis
        ]
        H_mpo = anyon_ham(model, sites)
        H_from_mpo = [
            real(inner(prime(bra), H_mpo, ket)) for bra in product_states,
            ket in product_states
        ]
        @test H_from_mpo ≈ anyon_ham(model) atol = 1e-13
    end
end

@testset "Measurement Operators" begin
    N = 4
    sites = siteinds("Qubit", N)
    τ = 1.0

    idx=2
    # Test measurement operator creation
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    M_p = measurement_operator_mps(model, sites, idx, τ, false)
    M_m = measurement_operator_mps(model, sites, idx, τ, true)

    s_im1 = sites[idx-1] # site 1
    s_i = sites[idx]   # site 2
    s_ip1 = sites[idx+1] # site 3
    row_inds = (prime(s_im1), prime(s_i), prime(s_ip1))
    col_inds = (s_im1, s_i, s_ip1)
    permuted_M_p = ITensors.permute(M_p, row_inds..., col_inds...)
    permuted_M_m = ITensors.permute(M_m, row_inds..., col_inds...)
    M_pmatrix = reshape(permuted_M_p.tensor.storage, 8, 8)
    M_mmatrix = reshape(permuted_M_m.tensor.storage, 8, 8)

    @test M_pmatrix^2 + M_mmatrix^2 ≈ I(8)
    # @test  M_pmatrix == 
    #       [1.0 0.0 0.0 0.0; 
    #        0.0 1.0 0.0 0.0; 
    #        0.0 0.0 exp(-τ) 0.0; 
    #        0.0 0.0 0.0 exp(-τ)]
    # @test M_mmatrix == 
    #       [exp(-τ) 0.0 0.0 0.0; 
    #        0.0 exp(-τ) 0.0 0.0; 
    #        0.0 0.0 1.0 0.0; 
    #        0.0 0.0 0.0 1.0]


    # Test invalid inputs
    @test_throws AssertionError measurement_operator_mps(model, sites, 0, τ, false)
    @test_throws AssertionError measurement_operator_mps(model, sites, N+1, τ, false)
end

@testset "Fibonacci PBC boundary MPO" begin
    N = 6
    τ = 1.0
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    ψ, sites = initial_mps(N)

    # The MPO boundary path must reproduce the former noncontiguous-gate path
    # without moving the boundary sites through the MPS with SWAP gates.
    for i in (1, N), sign in (false, true)
        gate = measurement_operator_mps(model, sites, i, τ, sign)
        ψ_gate, p_gate = FibonacciChain._measuremap_with_operator(
            ψ,
            gate;
            cutoff = 1e-14,
            maxdim = 1000,
        )
        ψ_mpo, p_mpo = measuremap(
            model,
            ψ,
            sites,
            i,
            τ,
            sign;
            cutoff = 1e-14,
            maxdim = 1000,
        )

        @test p_mpo ≈ p_gate atol = 1e-13
        @test abs(inner(ψ_mpo, ψ_gate))^2 ≈ 1 atol = 1e-13
    end
end

@testset "Fibonacci constraint projector MPO" begin
    N = 4
    sites = siteinds("Qubit", N)
    for pbc in (false, true)
        projector = FibonacciChain.fibonacci_constraint_projector_mpo(sites; pbc = pbc)
        for bits in Iterators.product(fill(0:1, N)...)
            state = productMPS(sites, collect(string.(bits)))
            projected = apply(projector, state; cutoff = 0.0)
            is_allowed = all(bits[i] == 0 || bits[i + 1] == 0 for i = 1:(N - 1))
            is_allowed &= !pbc || bits[N] == 0 || bits[1] == 0
            @test norm(projected) ≈ (is_allowed ? 1.0 : 0.0) atol = 1e-13
        end
    end
end

@testset "topological_charge_mpo" begin
    # The MPO must reproduce ⟨x|Y|x'⟩ = ∏ᵢ (F^{τ xᵢ τ}_{x'ᵢ₊₁})^{xᵢ₊₁}_{x'ᵢ}
    # for every pair of bit configurations, legal fusion path or not.
    state_mps(sites, N, buf) = productMPS(
        sites,
        [bit ? "1" : "0" for bit in reverse(digits(Bool, buf; base = 2, pad = N))],
    )

    for (N, pbc) in ((4, true), (4, false), (2, true))
        model = AnyonModel(FibonacciAnyon(), N; pbc = pbc)
        sites = siteinds("Qubit", N)
        Y_mpo = topological_charge_mpo(sites; pbc = pbc)
        for buf_in in 0:(2^N - 1)
            ψ_in = state_mps(sites, N, buf_in)
            mapped = apply(Y_mpo, ψ_in; cutoff = 1e-15, maxdim = 100)
            for buf_out in 0:(2^N - 1)
                ψ_out = state_mps(sites, N, buf_out)
                expected = Fsymmetry_coef(model, BitStr{N}(buf_out), BitStr{N}(buf_in))
                @test inner(ψ_out, mapped) ≈ expected atol = 1e-13
            end
        end
    end

    # On the constrained anyon basis the MPO must coincide with the exact
    # topological charge operator, whose PBC eigenvalues are ϕ (y=1 sector)
    # and -1/ϕ (y=τ sector).
    N = 8
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    Y_exact = topological_charge_operator(model)
    sites = siteinds("Qubit", N)
    Y_mpo = topological_charge_mpo(sites; pbc = true)
    basis = anyon_basis(model)
    basis_mps = [state_mps(sites, N, b.buf) for b in basis]
    Y_rebuilt = [
        inner(basis_mps[i], apply(Y_mpo, basis_mps[j]; cutoff = 1e-15, maxdim = 100))
        for i in eachindex(basis), j in eachindex(basis)
    ]
    @test Y_rebuilt ≈ Y_exact atol = 1e-10

    ϕ = (1 + √5) / 2
    y_eigs = eigvals((Y_rebuilt+Y_rebuilt')/2)
    @test all(y -> min(abs(y - ϕ), abs(y + inv(ϕ))) < 1e-8, y_eigs)
end

@testset "topological_charge_mpo bulk tensor and full contraction" begin
    ϕ = (1 + √5) / 2
    N = 4
    sites = siteinds("Qubit", N)

    # OBC bulk tensor: the automaton transition table
    # A[(xᵢ₋₁, x'ᵢ₋₁), (xᵢ, x'ᵢ) folded, (xᵢ, x'ᵢ)] = (F^{τ xᵢ₋₁ τ}_{x'ᵢ})^{xᵢ}_{x'ᵢ₋₁},
    # with the physical pair folded into the middle dim (prime index fastest).
    Y_obc = topological_charge_mpo(sites; pbc = false)
    A = reshape(
        Array(
            Y_obc[2],
            linkind(Y_obc, 1),
            prime(sites[2]),
            dag(sites[2]),
            linkind(Y_obc, 2),
        ),
        4, 4, 4,
    )
    expected = zeros(4, 4, 4)
    expected[1, 1, 1] = -ϕ^(-1)
    expected[2, 1, 1] = ϕ^(-1 / 2)
    expected[3, 1, 1] = 1.0
    expected[1, 3, 2] = 1.0
    expected[3, 3, 2] = 1.0
    expected[1, 2, 3] = ϕ^(-1 / 2)
    expected[2, 2, 3] = ϕ^(-1)
    @test A ≈ expected atol = 1e-14

    # Fully contracting the MPO must give ⟨x|Y|x'⟩ = Fsymmetry_coef(x, x')
    # for every pair of bit configurations, legal fusion path or not.
    pbc = true
    model = AnyonModel(FibonacciAnyon(), N; pbc = pbc); # obc has the same eigvals, but not the same matrix elements order
    Y_contracted = foldr(*, topological_charge_mpo(sites; pbc = pbc))
    i_in  = filter(i -> plev(i) > 0, inds(Y_contracted))
    i_out = filter(i -> plev(i) == 0, inds(Y_contracted))
    Y_perm = permute(Y_contracted, i_in..., i_out...)
    Y_matrix = reshape(Y_perm.tensor, 2^4, 2^4)
    idx = pbc ? [1,2,3,5,6,9,11] : [1,2,3,5,6,9,10,11]
    Y_ed = topological_charge_operator(model)
    @test Y_matrix[idx, idx] ≈ Y_ed    
end

@testset "topological_charge_mpo sector structure" begin
    ϕ = (1 + √5) / 2
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    sites = siteinds("Qubit", N)
    Y_mpo = topological_charge_mpo(sites; pbc = true)
    basis = anyon_basis(model)

    basis_product_mps(buf) = productMPS(
        sites,
        [bit ? "1" : "0" for bit in reverse(digits(Bool, buf; base = 2, pad = N))],
    )
    # Random states supported on the constrained anyon basis
    rng = MersenneTwister(42)
    function random_physical_state()
        coef = randn(rng, length(basis))
        ψ = coef[1] * basis_product_mps(basis[1].buf)
        for i in 2:length(basis)
            ψ = add(ψ, coef[i] * basis_product_mps(basis[i].buf); cutoff = 1e-15, maxdim = 200)
        end
        normalize!(ψ)
        return ψ
    end
    # Sector projectors via the Y MPO: P1 = (Y + ϕ⁻¹ I)/(ϕ + ϕ⁻¹),
    # Pτ = (Y - ϕ I)/(-ϕ⁻¹ - ϕ)
    project_1(ψ) = inv(ϕ + inv(ϕ)) * add(
        apply(Y_mpo, ψ; cutoff = 1e-15, maxdim = 200),
        inv(ϕ) * ψ;
        cutoff = 1e-15,
        maxdim = 200,
    )
    project_τ(ψ) = inv(-inv(ϕ) - ϕ) * add(
        apply(Y_mpo, ψ; cutoff = 1e-15, maxdim = 200),
        (-ϕ) * ψ;
        cutoff = 1e-15,
        maxdim = 200,
    )

    ψ = random_physical_state()
    φ = random_physical_state()
    Yψ = apply(Y_mpo, ψ; cutoff = 1e-15, maxdim = 200)

    # On the constrained subspace Y has eigenvalues ϕ and -1/ϕ, hence the
    # minimal polynomial Y² = Y + I
    YYψ = apply(Y_mpo, Yψ; cutoff = 1e-15, maxdim = 200)
    residual = add(
        YYψ,
        (-1.0) * add(Yψ, ψ; cutoff = 1e-15, maxdim = 200);
        cutoff = 1e-15,
        maxdim = 200,
    )
    @test norm(residual) < 1e-10

    # Y is hermitian: ⟨ψ|Y|φ⟩ = ⟨φ|Y|ψ⟩*
    @test inner(prime(ψ), Y_mpo, φ) ≈ conj(inner(prime(φ), Y_mpo, ψ)) atol = 1e-12

    # The sector projectors are complete and orthogonal on physical states
    p1ψ = project_1(ψ)
    pτψ = project_τ(ψ)
    @test norm(add(p1ψ, pτψ; cutoff = 1e-15, maxdim = 200) - ψ) < 1e-10
    @test norm(project_1(pτψ)) < 1e-10
    @test real(inner(prime(p1ψ), Y_mpo, p1ψ)) / inner(p1ψ, p1ψ) ≈ ϕ atol = 1e-10
    @test real(inner(prime(pτψ), Y_mpo, pτψ)) / inner(pτψ, pτψ) ≈ -inv(ϕ) atol = 1e-10

    # Y commutes with every measurement transfer matrix: [T, Y] = 0 for a
    # fixed-outcome period applied with normalized = false
    τ = 0.8
    n_layers = FibonacciChain.layers_per_period(model)
    sample = BitMatrix(ones(Int8, n_layers, FibonacciChain._samples_per_layer(model)))
    sample[1, 1] = 0
    config = MeasureConfig(τ = τ, mode = :sample, t₂ = 1, enable_τ_eff = false)
    apply_T(state) = FibonacciChain._sample_measure_mps(
        model,
        sites,
        state,
        sample,
        config;
        cutoff = 1e-15,
        maxdim = 200,
        normalized = false,
    ).state
    TYψ = apply_T(Yψ)
    YTψ = apply(Y_mpo, apply_T(ψ); cutoff = 1e-15, maxdim = 200)
    @test norm(add(TYψ, (-1.0) * YTψ; cutoff = 1e-15, maxdim = 200)) / norm(TYψ) < 1e-10

    # The automaton bond dimensions: 4 for open boundaries (carrying the
    # current fusion pair), 16 for periodic boundaries (also closing the ring)
    Y_obc = topological_charge_mpo(sites; pbc = false)
    @test maxlinkdim(Y_obc) == 4
    @test maxlinkdim(Y_mpo) == 16
end

@testset "TCI ground state is a y=1 eigenstate of the Y MPO" begin
    # The periodic Fibonacci (tricritical Ising) ground state carries trivial
    # topological charge, so Y|GS⟩ = ϕ|GS⟩. The exact ground state is used
    # because penalty-based DMRG mixes the nearly-degenerate sectors at the
    # 1e-3 level; the state is assembled as an MPS and probed purely through
    # the Y MPO.
    ϕ = (1 + √5) / 2
    L = 10
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    basis = anyon_basis(model)
    gs_vec = real(Arpack.eigs(anyon_ham_sparse(model); nev = 1, which = :SR)[2][:, 1])

    sites = siteinds("Qubit", L)
    product_state(i) = productMPS(
        sites,
        [bit ? "1" : "0" for bit in reverse(digits(Bool, basis[i].buf; base = 2, pad = L))],
    )
    gs = gs_vec[1] * product_state(1)
    for i in 2:length(basis)
        gs = add(gs, gs_vec[i] * product_state(i); cutoff = 1e-15, maxdim = 256)
    end
    normalize!(gs)

    Y_mpo = topological_charge_mpo(sites; pbc = true)
    y_exp = real(inner(prime(gs), Y_mpo, gs))
    @test y_exp ≈ ϕ atol = 1e-8

    Yψ = apply(Y_mpo, gs; cutoff = 1e-15, maxdim = 256)
    @test norm(add(Yψ, (-ϕ) * gs; cutoff = 1e-15, maxdim = 512)) < 1e-7

    # The penalty DMRG ground state at L = 10 (20 sweeps, maxdim 128,
    # seed 1234) is already a y=1 eigenstate at the 1e-3 level, unlike the
    # nearly-degenerate L = 16 case where sector mixing is visible.
    # Note the DMRG MPS carries its own site indices, so the Y MPO must be
    # rebuilt on `siteinds(ψ)`.
    ψ, _ = anyon_mps_gst(
        model;
        sweep_times = 20,
        maxdim = 128,
        cutoff = 1e-12,
        outputlevel = 0,
        seed = 1234,
    )
    Y_mpo_dmrg = topological_charge_mpo(siteinds(ψ); pbc = true)
    y_dmrg = real(inner(prime(ψ), Y_mpo_dmrg, ψ)) / real(inner(ψ, ψ))
    @test y_dmrg ≈ ϕ atol = 1e-5  # measured: 1.61803253

    Yψ_dmrg = apply(Y_mpo_dmrg, ψ; cutoff = 1e-15, maxdim = 512)
    res_dmrg = norm(add(Yψ_dmrg, (-ϕ) * ψ; cutoff = 1e-15, maxdim = 1024))
    @test res_dmrg < 3e-3  # measured: 1.8e-3
end

@testset "Potts ground state is a y=τ eigenstate of the Y MPO" begin
    # The ferromagnetically coupled periodic Fibonacci chain (3-state Potts
    # critical point) has a unique ground state carrying total charge τ,
    # so Y|GS⟩ = -ϕ⁻¹|GS⟩. Exact ground state assembled as an MPS, probed
    # purely through the Y MPO (mirrors the TCI testset above).
    ϕ = (1 + √5) / 2
    L = 10
    model = AnyonModel(FibonacciAnyon(), L; pbc = true, measure_operator = :Ferro)
    basis = anyon_basis(model)
    gs_vec = real(Arpack.eigs(anyon_ham_sparse(model); nev = 1, which = :SR)[2][:, 1])

    sites = siteinds("Qubit", L)
    product_state(i) = productMPS(
        sites,
        [bit ? "1" : "0" for bit in reverse(digits(Bool, basis[i].buf; base = 2, pad = L))],
    )
    gs = gs_vec[1] * product_state(1)
    for i in 2:length(basis)
        gs = add(gs, gs_vec[i] * product_state(i); cutoff = 1e-15, maxdim = 256)
    end
    normalize!(gs)

    Y_mpo = topological_charge_mpo(sites; pbc = true)
    y_exp = real(inner(prime(gs), Y_mpo, gs))
    @test y_exp ≈ -inv(ϕ) atol = 1e-8

    Yψ = apply(Y_mpo, gs; cutoff = 1e-15, maxdim = 256)
    @test norm(add(Yψ, inv(ϕ) * gs; cutoff = 1e-15, maxdim = 512)) < 1e-7

    # Penalty DMRG ground state at L = 10 (20 sweeps, maxdim 128, seed 1234):
    # a y=τ eigenstate at the 1e-3 level. The DMRG MPS carries its own site
    # indices, so the Y MPO must be rebuilt on `siteinds(ψ)`.
    ψ, _ = anyon_mps_gst(
        model;
        sweep_times = 20,
        maxdim = 128,
        cutoff = 1e-12,
        outputlevel = 0,
        seed = 1234,
    )
    Y_mpo_dmrg = topological_charge_mpo(siteinds(ψ); pbc = true)
    y_dmrg = real(inner(prime(ψ), Y_mpo_dmrg, ψ)) / real(inner(ψ, ψ))
    @test y_dmrg ≈ -inv(ϕ) atol = 1e-5  # measured: -0.61803281

    Yψ_dmrg = apply(Y_mpo_dmrg, ψ; cutoff = 1e-15, maxdim = 512)
    res_dmrg = norm(add(Yψ_dmrg, inv(ϕ) * ψ; cutoff = 1e-15, maxdim = 1024))
    @test res_dmrg < 3e-3  # measured: 1.6e-3
end

@testset "Constraint-preserving post-selection MPS" begin
    N = 8
    periods = 2N
    τ = atanh(0.95)
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    samples = BitMatrix(ones(Bool, 2periods, N ÷ 2))

    exact_state = zeros(length(anyon_basis(model)))
    exact_state[1] = 1.0
    exact = bulk_evolution(
        model,
        exact_state,
        MeasureConfig(τ = τ, mode = :sample, t₂ = periods),
        samples,
    )

    ψ, sites = initial_mps(N)
    mps = bulk_evolution(
        model,
        sites,
        ψ,
        MeasureConfig(
            τ = τ,
            mode = :sample,
            t₂ = periods,
            cutoff = 1e-12,
            maxdim = 100,
            enforce_fibonacci_constraint = true,
        ),
        samples,
    )

    @test mps.entanglement_entropys ≈ exact.entanglement_entropys atol = 2e-6
    @test anyon_eelis(model, mps.state) ≈ anyon_eelis(model, exact.state) atol = 2e-6
end

@testset "Single Measurement Application" begin
    N = 4
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    τ = 1.0

    sites = siteinds("Qubit", N)

    # Create initial product state (vacuum state)
    state = ["0" for _ = 1:N]

    # Create MPS from product state
    ψ0 = random_mps(sites, state)

    # Apply measurements
    ψ_p, prob_p = measuremap(model, ψ0, sites, 1, τ, false;)
    ψ_m, prob_m = measuremap(model, ψ0, sites, 1, τ, true;)

    st = zeros(length(anyon_basis(model)));
    st[1] = 1.0
    state_after_p = measuremap(model, τ, st, 1, false)
    p = state_after_p'*state_after_p

    @test prob_p ≈ p
    @test prob_m ≈ 1 - p  # Should be orthogonal to ψ_p

    if prob_p > 1e-12
        @test abs(inner(ψ_p, ψ_p) - 1.0) < 1e-10
    end
    if prob_m > 1e-12
        @test abs(inner(ψ_m, ψ_m) - 1.0) < 1e-10
    end
end


@testset "Measurement Enumeration" begin
    N = 4
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    τ = 1.0

    sites = siteinds("Qubit", N)

    # Create initial product state (vacuum state)
    state = ["0" for _ = 1:N]
    state_exact = zeros(length(anyon_basis(model)))
    state_exact[1] = 1.0  # Vacuum state
    # Use minimal measurement sites for enumeration
    measurement_sites = collect(2:2:N)

    ψ = random_mps(sites, state)
    # Enumerate trajectories
    final_states, trajectories, probabilities =
        mps_measurement_enumeration(model, ψ, sites, measurement_sites, τ;)

    # final_states_vector = map(x-> reduce(*, x).tensor.storage, final_states) # Noting its elements arranging is different from the final_states_exact, they choose different definition?

    final_states_exact, trajectories_exact, probabilities_exact =
        measurement_enumeration(model, τ, state_exact, measurement_sites)

    @test trajectories == trajectories_exact
    @test probabilities ≈ probabilities_exact

    # Check probability normalization
    @test abs(sum(probabilities) - 1.0) < 1e-10

    # Check that all probabilities are positive
    @test all(p -> p > 0, probabilities)
end

@testset "Boundary_Born" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    τ = 1.0

    ψ, sites = initial_mps(N)

    seed=10
    # Perform sampling
    st = zeros(length(anyon_basis(model)));
    st[1] = 1.0

    measure_config = MeasureConfig(
        τ = τ,
        t₂ = 1,
        rng = MersenneTwister(seed),
        mode = :Born,
        enable_τ_eff = false,
    )
    measure_outcome_mps = boundary_evolution(model, sites, ψ, measure_config)
    sample_measured_states_mps, sample_mps, samples_free_energy_mps =
        measure_outcome_mps.state,
        measure_outcome_mps.sample,
        measure_outcome_mps.free_energy

    measure_config = MeasureConfig(
        τ = τ,
        t₂ = 1,
        rng = MersenneTwister(seed),
        mode = :Born,
        enable_τ_eff = false,
    ) # NEED to reset rng to ensure same sampling
    measure_outcome = boundary_evolution(model, st, measure_config)
    sample_measured_states, sample, samples_free_energy =
        measure_outcome.state, measure_outcome.sample, measure_outcome.free_energy

    @test sample_mps == sample
    @test samples_free_energy_mps ≈ samples_free_energy
end

@testset "_born_measure" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N)
    t = 10
    measure_config =
        MeasureConfig(τ = 1000.0, t₂ = t, rng = MersenneTwister(42), mode = :Born)

    state = zeros(length(anyon_basis(model)));
    state[1] = 1.0
    measure_outcome = FibonacciChain._born_measure(model, state, measure_config)
    samples, sample_free_energy = measure_outcome.samples, measure_outcome.free_energys
    @test size(samples) == (20, 3)
    @test sample_free_energy[end] ≈ 1.5009765892377303 atol=1e-6

    measure_config =
        MeasureConfig(τ = 1000.0, t₂ = t, rng = MersenneTwister(42), mode = :Born)
    ψ, sites = initial_mps(N)
    measure_outcome_mps = FibonacciChain._born_measure_mps(model, sites, ψ, measure_config)
    samples_mps, sample_free_energy_mps =
        measure_outcome_mps.samples, measure_outcome_mps.free_energys
    @test samples_mps == samples
    @test sample_free_energy_mps ≈ sample_free_energy

    measure_config2 = MeasureConfig(
        τ = 1000.0,
        t₂ = t,
        rng = MersenneTwister(42),
        mode = :Born,
        truncate_every_events = 2,
    )
    measure_outcome_mps2 =
        FibonacciChain._born_measure_mps(model, sites, ψ, measure_config2)
    samples_mps2, sample_free_energy_mps2 =
        measure_outcome_mps2.samples, measure_outcome_mps2.free_energys
    @test samples_mps2 == samples
    @test sample_free_energy_mps2 ≈ sample_free_energy
end

@testset "_sample_measure" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N)
    t = 10
    measure_config = MeasureConfig(τ = 1000.0, t₂ = t, mode = :sample)
    state = zeros(length(anyon_basis(model)));
    state[1] = 1.0
    samples = BitMatrix(undef, 2t, div(N, 2))

    measure_outcome = FibonacciChain._sample_measure(model, state, samples, measure_config)
    samples, sample_free_energy = measure_outcome.samples, measure_outcome.free_energys
    @test size(samples) == (20, 3)
    @test sample_free_energy[end] ≈ 0.5385529416309107 atol=1e-6

    ψ, sites = initial_mps(N)
    measure_outcome_mps =
        FibonacciChain._sample_measure_mps(model, sites, ψ, samples, measure_config)
    samples_mps, sample_free_energy_mps =
        measure_outcome_mps.samples, measure_outcome_mps.free_energys
    @test samples_mps == samples
    @test sample_free_energy_mps ≈ sample_free_energy
end

@testset "Bulk_Born" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    τ = 1.0
    D = 2  # Small number of layers for testing 

    ψ, sites = initial_mps(N)

    seed=10
    # Perform sampling
    st = zeros(length(anyon_basis(model)));
    st[1] = 1.0

    measure_config = MeasureConfig(τ = τ, t₂ = D, rng = MersenneTwister(seed), mode = :Born)
    # Perform bulk measurements
    measure_outcome_mps = bulk_evolution(model, sites, ψ, measure_config)
    bulk_samples, bulk_free_energy =
        measure_outcome_mps.samples, measure_outcome_mps.free_energys

    measure_config = MeasureConfig(τ = τ, t₂ = D, rng = MersenneTwister(seed), mode = :Born) # NEED to reset rng to ensure same sampling
    measure_outcome = bulk_evolution(model, st, measure_config)
    bulk_samples_exact, bulk_free_energy_exact =
        measure_outcome.samples, measure_outcome.free_energys

    @test bulk_samples == bulk_samples_exact
    @test bulk_free_energy ≈ bulk_free_energy_exact
end

@testset "_apply_measurement_layer" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)
    τ = 1.0

    ψ, sites = initial_mps(N)
    st = zeros(length(anyon_basis(model)));
    st[1] = 1.0

    # Apply measurement to a specific layer
    measurement_layer = 2
    bulk_samples = BitVector(ones(3))

    measure_outcome_mps = FibonacciChain._apply_measurement_layer_mps(
        model,
        τ,
        sites,
        ψ,
        bulk_samples,
        measurement_layer,
    )
    ψ_layer, F = measure_outcome_mps.state, measure_outcome_mps.free_energy

    measure_outcome = FibonacciChain._apply_measurement_layer(
        model,
        τ,
        st,
        bulk_samples,
        layer_idx = measurement_layer,
    )
    st_exact, F_exact = measure_outcome.state, measure_outcome.free_energy
    # Here we use keyword argument to avoid confusion

    inds = [i.buf for i in anyon_basis(model)] .+ 1

    ψ_dense = reduce(*, ψ_layer).tensor.storage[inds]
    @test ψ_dense[ψ_dense .> 0] ≈ st_exact[st_exact .> 0]
    @test F ≈ F_exact
end

@testset "bulk_evolution" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)

    ψ, sites = initial_mps(N)
    st = zeros(length(anyon_basis(model)));
    st[1] = 1.0

    # Generate a specific state

    τ = 1.0  # Example τ value
    bulk_samples = BitMatrix([1 1 1; 0 0 0])
    measure_config = MeasureConfig(τ = τ, t₂ = 1, mode = :sample)
    measure_outcome_mps = bulk_evolution(model, sites, ψ, measure_config, bulk_samples)
    generated_statelis, F = measure_outcome_mps.state, measure_outcome_mps.free_energys
    measure_outcome = bulk_evolution(model, st, measure_config, bulk_samples)
    generated_statelis_exact, F_exact = measure_outcome.state, measure_outcome.free_energys

    inds = [i.buf for i in anyon_basis(model)] .+ 1

    # Convert MPS states to dense vectors for comparison, note the order of elements may differ, as actually we are dealing with OBC MPS, but PBC Hamiltonian.
    ψ_dense = reduce(*, generated_statelis).tensor.storage[inds]
    ψ_dense = sort(ψ_dense[ψ_dense .> 0])
    generated_statelis_exact = sort(generated_statelis_exact[generated_statelis_exact .> 0])

    @test ψ_dense ≈ generated_statelis_exact
    @test F ≈ F_exact
end

@testset "Entanglement Entropy" begin
    N = 6
    model = AnyonModel(FibonacciAnyon(), N; pbc = true)

    sites = siteinds("Qubit", N)

    # Create initial product state (vacuum state)
    state = ["0" for _ = 1:N]

    ψ = random_mps(sites, state)
    # Calculate entanglement entropy at different cuts
    EElis = anyon_eelis(model, ψ)
    @test all(EElis .>= 0)  # Entanglement entropy should be non-negative
end

function samples_generate_mps(L::Int64, τ::Float64, seed::Int64, D::Int64 = 5L)
    rng = MersenneTwister(seed)

    ψ, sites = initial_mps(L)

    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    measure_config = MeasureConfig(τ = τ, t₂ = D, rng = rng, mode = :Born)
    measure_outcome = bulk_evolution(model, sites, ψ, measure_config)
    sample, sample_free_energy = measure_outcome.samples, measure_outcome.free_energys

    halfchain_EE_tlis = measure_outcome.entanglement_entropys
    final_state = measure_outcome.state
    final_EElis = anyon_eelis(model, final_state)

    return sample, sample_free_energy, final_EElis, halfchain_EE_tlis
end

function samples_generate(L::Int64, τ::Float64, seed::Int64, D::Int64 = 5L)
    rng = MersenneTwister(seed)

    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    measure_config = MeasureConfig(τ = τ, t₂ = D, rng = rng, mode = :Born)
    st = zeros(length(anyon_basis(model)))
    st[1] = 1.0

    measure_outcome = bulk_evolution(model, st, measure_config)
    sample, sample_free_energy = measure_outcome.samples, measure_outcome.free_energys

    halfchain_EE_tlis = measure_outcome.entanglement_entropys
    final_state = measure_outcome.state
    final_EElis = anyon_eelis(model, final_state)

    return sample, sample_free_energy, final_EElis, halfchain_EE_tlis
end

@testset "Observable" begin
    L = 6
    τ = 1.0
    D = 5

    # Generate samples
    seed = 42
    sample, sample_free_energy, final_EElis, halfchain_EE_tlis =
        samples_generate(L, τ, seed, D)

    sample_mps, sample_free_energy_mps, final_EElis_mps, halfchain_EE_tlis_mps =
        samples_generate_mps(L, τ, seed, D)

    @test sample == sample_mps
    @test sample_free_energy ≈ sample_free_energy_mps
    @test final_EElis ≈ final_EElis_mps
    @test halfchain_EE_tlis ≈ halfchain_EE_tlis_mps
end

@testset "lyapunov_spectrum_mps" begin
    L = 8
    τ = atanh(0.95)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    sites = siteinds("Qubit", L)

    # Reference: direct diagonalization of a single-period transfer matrix
    sample_ref = BitMatrix(ones(Int8, 2, div(L, 2)))
    T = transfer_matrix(model, τ, sample_ref)
    energy = eigen(T).values
    sorted_energy = sort(energy, by = abs, rev = true)
    spectrum_ref = -log.(real.(sorted_energy[1:10]))

    # Multi-period sample: subspace iteration should converge for dominant eigenvalues
    D = 50
    sample_long = BitMatrix(ones(Int8, D, div(L, 2)))
    spectrum_sub = lyapunov_spectrum_mps(
        model, sites, τ, sample_long; n_states = 10, cutoff = 1e-12, maxdim = 100
    )
    @test size(spectrum_sub) == (10, div(D, 2))
    # The dominant eigenvalues converge well; subdominant ones are more sensitive
    # to MPS truncation and orthogonalization errors.
    @test spectrum_sub[1:3, end] ≈ spectrum_ref[1:3] atol = 1e-5
end

@testset "lyapunov_spectrum_mps" begin
    L = 8
    τ = atanh(0.95)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    sites = siteinds("Qubit", L)

    # Reference: direct diagonalization of a single-period transfer matrix
    D = 20
    binary_distribution(p, rng) = rand(rng) < p ? 1 : 0
    rng = MersenneTwister(42)
    sample_ref = BitMatrix(
                reshape([binary_distribution(0.4, rng) for _ = 1:D*div(L, 2)], D, div(L, 2)))
    spectrum_ref = lyapunov_spectrum(model, τ, sample_ref)

    # Multi-period sample: subspace iteration should converge for dominant eigenvalues
    rng = MersenneTwister(42)
    sample_long = BitMatrix(
                reshape([binary_distribution(0.4, rng) for _ = 1:D*div(L, 2)], D, div(L, 2)))
    spectrum_sub = lyapunov_spectrum_mps(
        model, sites, τ, sample_long; n_states = 10, cutoff = 1e-12, maxdim = 100
    )
    @test size(spectrum_sub) == (10, div(D, 2))
    # The dominant eigenvalues converge well; subdominant ones are more sensitive
    # to MPS truncation and orthogonalization errors.
    @test spectrum_sub[1:7, :] ≈ spectrum_ref[1:7, :] atol = 1e-5
end


@testset "Y MPO sector tools: initial state, tracking, sector Lyapunov" begin
    L = 8
    τ = atanh(0.95)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    sites = siteinds("Qubit", L)
    ϕ = (1 + √5) / 2
    Y_mpo = topological_charge_mpo(sites; pbc = true)

    basis = anyon_basis(model)
    basis_mps(buf) = productMPS(
        sites,
        [bit ? "1" : "0" for bit in reverse(digits(Bool, buf; base = 2, pad = L))],
    )

    # y=1 projected all-τ state built with the Y MPO: P1|0⋯0⟩ ∝ (Y + ϕ⁻¹ I)|0⋯0⟩
    vacuum = basis_mps(0)
    y1_init = add(
        apply(Y_mpo, vacuum; cutoff = 1e-14, maxdim = 100),
        inv(ϕ) * vacuum;
        cutoff = 1e-14,
        maxdim = 100,
    )
    normalize!(y1_init)

    # Compare against the exact sector projection
    Y_exact = topological_charge_operator(model)
    P1 = (Y_exact + inv(ϕ) * I) / (ϕ + inv(ϕ))
    exact_init = zeros(length(basis))
    exact_init[1] = 1.0
    exact_init = P1 * exact_init
    exact_init ./= norm(exact_init)
    mps_amplitudes = [inner(basis_mps(b.buf), y1_init) for b in basis]
    @test mps_amplitudes ≈ exact_init atol = 1e-10

    # ⟨Y⟩ via the MPO equals the y=1 eigenvalue ϕ
    y_exp =
        real(inner(prime(y1_init), Y_mpo, y1_init)) / real(inner(y1_init, y1_init))
    @test y_exp ≈ ϕ atol = 1e-10

    # Born evolution tracks the y=1 sector on-the-fly via the Y MPO
    periods = 5
    outcome = bulk_evolution(
        model,
        sites,
        y1_init,
        MeasureConfig(
            τ = τ,
            t₂ = periods,
            mode = :Born,
            rng = MersenneTwister(1),
            enable_τ_eff = false,
            cutoff = 1e-12,
            maxdim = 100,
            track_y_expectation = true,
        ),
    )
    @test length(outcome.y_expectation_values) == periods
    # y_expectation_values are stored as Float32, so the tolerance must exceed
    # Float32 roundoff (~1.2e-7)
    @test all(y -> abs(y - ϕ) < 1e-6, outcome.y_expectation_values)

    # Sector-restricted Lyapunov spectrum with an MPS frame: project the first
    # basis product states with P1 via the Y MPO, then let `initial_states`
    # drive the subspace iteration. `sector = :trivial` re-projects the frame
    # into y=1 after every period; without it, MPS truncation noise leaks into
    # the y=τ sector, whose leading exponent (-2.95 here) outgrows the y=1
    # subleading directions. Compare at the same time horizon against
    # the exact sector Lyapunov function, so only MPS truncation error remains.
    n_states = 3
    frame = Vector{MPS}(undef, n_states)
    for i in 1:n_states
        product_state = basis_mps(basis[i].buf)
        frame[i] = add(
            apply(Y_mpo, product_state; cutoff = 1e-12, maxdim = 100),
            inv(ϕ) * product_state;
            cutoff = 1e-12,
            maxdim = 100,
        )
    end

    D = 400
    sample_long = BitMatrix(ones(Int8, D, div(L, 2)))
    spectrum_mps = lyapunov_spectrum_mps(
        model,
        sites,
        τ,
        sample_long;
        initial_states = frame,
        sector = :trivial,
        cutoff = 1e-12,
        maxdim = 100,
    )
    lyapunov_mps =
        cumsum(-spectrum_mps; dims = 2) ./ reshape(collect(1:div(D, 2)), 1, :)

    exact = lyapunov_spectrum_topological_sector(
        model,
        τ,
        sample_long;
        sector = :trivial,
        n_states = n_states,
    )
    @test lyapunov_mps[:, end] ≈ exact.lyapunov_exponents[:, end] atol = 1e-3
end
