using FibonacciChain, Test, LinearAlgebra, SparseArrays, Random

@testset "Periodic Ising fusion Hilbert space" begin
    for L in 1:10
        model = AnyonModel(IsingAnyon(), L)
        basis = anyon_basis(model)
        @test length(basis) == (isodd(L) ? 0 : 2^(L÷2+1))
        @test issorted(basis) && allunique(basis)
        @test all(s -> all(i -> xor(s[i] == 2, s[mod1(i+1,L)] == 2), 1:L), basis)
        if iseven(L)
            @test count(s -> s[1] == 2, basis) == 2^(L÷2)
        else
            @test size(anyon_ham(model)) == (0, 0)
        end
    end
    @test_throws ArgumentError AnyonModel(IsingAnyon(), 0)
    @test_throws ArgumentError AnyonModel(IsingAnyon(), 4; pbc=false)
    @test_throws ArgumentError AnyonModel(IsingAnyon(), 4; measure_operator=:X)
    @test_throws ArgumentError AnyonModel(IsingAnyon(), 4; J=Inf)
    @test_throws ArgumentError anyon_basis(AnyonModel(IsingAnyon(), 4); symmetry_block=:I)
end

@testset "Ising F move, projectors and Hamiltonian" begin
    F = [1 1; 1 -1]/sqrt(2)
    model = AnyonModel(IsingAnyon(), 6)
    basis = anyon_basis(model)
    # Direct F† |I><I| F comparison on a fixed σ–(I/η)–σ neighborhood.
    inds = [findfirst(==(Tuple(UInt8.(s))), basis) for s in
        ((2,0,2,0,2,0), (2,1,2,0,2,0))]
    @test Matrix(ising_fusion_projector(model, 2))[inds, inds] ≈ F'*[1 0; 0 0]*F
    for L in (2, 4, 6, 8)
        model = AnyonModel(IsingAnyon(), L; J=1.3, h=0.7)
        basis = anyon_basis(model)
        Ps = [ising_fusion_projector(model, i) for i in 1:L]
        for i in 1:L
            P, Q = Ps[i], ising_fusion_projector(model, i; channel=:eta)
            @test ishermitian(P) && ishermitian(Q)
            @test P*P ≈ P
            @test Q*Q ≈ Q
            @test P+Q ≈ I
            @test norm(P*Q) < 1e-12
        end
        H = anyon_ham(model)
        @test ishermitian(H)
        @test H ≈ Matrix(anyon_ham_sparse(model))
        @test H ≈ -sum((isodd(i) ? 1.3 : 0.7)*Ps[i] for i in 1:L)
        # Each sector is an independently constructed TFIM with the exact
        # projector shift and factor; staggered couplings swap between sectors.
        for parity in (false, true)
            sector = findall(s -> (s[1] == 2) == parity, basis)
            spin = AnyonModel(SpinHalf(), L÷2; model_type=:Ising,
                J=parity ? 1.3 : 0.7, h=parity ? 0.7 : 1.3)
            @test H[sector,sector] ≈ anyon_ham(spin)/2 - (L÷2)*(1.3+0.7)/2*I
            @test iszero(H[sector,setdiff(eachindex(basis),sector)])
        end
    end
    @test_throws ArgumentError ising_fusion_projector(model, 0)
    @test_throws ArgumentError ising_fusion_projector(model, 1; channel=:sigma)
end

@testset "Periodic Ising Temperley–Lieb algebra" begin
    # For either fixed fusion channel, e_i = √2 P_i and loop weight d = √2.
    # L=2 identifies the two neighbors and is not a distinct three-anyon test.
    for L in (4, 6, 8, 10), channel in (:I, :eta)
        @testset "L=$L, channel=$channel" begin
            model = AnyonModel(IsingAnyon(), L)
            Ps = [ising_fusion_projector(model, i; channel=channel) for i in 1:L]
            es = sqrt(2) .* Ps
            for i in 1:L
                @test Ps[i]^2 ≈ Ps[i] atol=1e-12 rtol=0
                @test es[i]^2 ≈ sqrt(2)*es[i] atol=1e-12 rtol=0
                for j in (mod1(i-1, L), mod1(i+1, L))
                    @test Ps[i]*Ps[j]*Ps[i] ≈ Ps[i]/2 atol=1e-12 rtol=0
                    @test es[i]*es[j]*es[i] ≈ es[i] atol=1e-12 rtol=0
                end
                for j in (i+1):L
                    distance = min(j-i, L-(j-i))
                    if distance > 1
                        @test Ps[i]*Ps[j] ≈ Ps[j]*Ps[i] atol=1e-12 rtol=0
                        @test es[i]*es[j] ≈ es[j]*es[i] atol=1e-12 rtol=0
                    end
                end
            end
        end
    end
end

@testset "Ising fusion projectors commute with topological loops" begin
    # Independent loop construction from the Ising TY+ F-symbol tensor:
    # <x'|Yσ|x> = ∏_i T[(x_i,x'_i),(x_{i+1},x'_{i+1})].
    # Pair order: (I,σ), (σ,I), (σ,ψ), (ψ,σ); code labels I=0, ψ=1, σ=2.
    # This oracle uses neither the local projector implementation nor a
    # single-spin-lattice D (which is not the full Hermitian Yσ).
    pair_index = Dict((0,2)=>1, (2,0)=>2, (2,1)=>3, (1,2)=>4)
    Tσ = [0.0 1.0 1.0 0.0;
          inv(sqrt(2)) 0.0 0.0 inv(sqrt(2));
          inv(sqrt(2)) 0.0 0.0 -inv(sqrt(2));
          0.0 1.0 -1.0 0.0]
    for L in (2, 4, 6, 8, 10)
        @testset "L=$L" begin
            model = AnyonModel(IsingAnyon(), L)
            basis = anyon_basis(model)
            dim = length(basis)
            Yσ, Yψ = zeros(dim, dim), zeros(dim, dim)
            for (col, ket) in enumerate(basis)
                flipped = map(x -> x == 2 ? x : UInt8(1)-x, ket)
                Yψ[searchsortedfirst(basis, flipped), col] = 1
                for (row, bra) in enumerate(basis)
                    pairs = [get(pair_index, (ket[i], bra[i]), 0) for i in 1:L]
                    if all(!iszero, pairs)
                        Yσ[row, col] = prod(Tσ[pairs[i], pairs[mod1(i+1,L)]] for i in 1:L)
                    end
                end
            end

            # Validate the oracle, including nonzero rank and sector exchange,
            # so a vanishing or identity loop cannot pass the commutator tests.
            Id = Matrix{Float64}(I, dim, dim)
            @test Yσ ≈ Yσ' atol=1e-12 rtol=0
            @test Yψ ≈ Yψ' atol=1e-12 rtol=0
            @test Yψ^2 ≈ Id atol=1e-12 rtol=0
            @test Yσ^2 ≈ Id+Yψ atol=1e-12 rtol=0
            @test Yσ*Yψ ≈ Yσ atol=1e-12 rtol=0
            @test Yψ*Yσ ≈ Yσ atol=1e-12 rtol=0
            @test rank(Yσ; atol=1e-12) == 2^(L÷2)
            @test tr(Yσ) ≈ 0 atol=1e-12
            @test tr(Yψ) ≈ 0 atol=1e-12
            A = findall(s -> s[1] == 2, basis)
            B = findall(s -> s[1] != 2, basis)
            @test iszero(Yσ[A,A]) && iszero(Yσ[B,B])
            @test iszero(Yψ[A,B]) && iszero(Yψ[B,A])

            # Both local fusion channels, at every link including the seam.
            for i in 1:L, channel in (:I, :eta)
                P = ising_fusion_projector(model, i; channel=channel)
                @test P*Yσ ≈ Yσ*P atol=1e-12 rtol=0
                @test P*Yψ ≈ Yψ*P atol=1e-12 rtol=0
            end
        end
    end
end

@testset "Ising fusion Kraus measurements" begin
    model = AnyonModel(IsingAnyon(), 6)
    basis = anyon_basis(model)
    rng = MersenneTwister(7)
    state = normalize(randn(rng, ComplexF64, length(basis)))
    for i in 1:6, τ in (0.0, 0.6, 1000.0, Inf)
        M0 = FibonacciChain.measure_matrix(model, τ, i, false)
        M1 = FibonacciChain.measure_matrix(model, τ, i, true)
        @test M0'*M0 + M1'*M1 ≈ I
        @test norm(M0*state)^2 + norm(M1*state)^2 ≈ 1
        for sign in (false,true)
            M = sign ? M1 : M0
            @test measuremap(model, τ, state, i, sign) ≈ M*state
            for (col,s) in enumerate(basis)
                r = measure_basismap(model, τ, s, i, sign)
                v = zeros(length(basis))
                v[searchsortedfirst(basis,r.s1)] += r.w1
                v[searchsortedfirst(basis,r.s2)] += r.w2
                @test v ≈ M[:,col]
            end
        end
        P = Matrix(ising_fusion_projector(model,i))
        if τ == 0
            @test M0 ≈ Matrix{Float64}(I,length(basis),length(basis))/sqrt(2)
        elseif isinf(τ)
            @test M0 ≈ P
            @test M1 ≈ I-P
        elseif τ < 1
            @test M0 ≈ exp(τ/2*(2P-I))/sqrt(2cosh(τ))
        end
    end
    @test_throws ArgumentError FibonacciChain.measure_matrix(model, -1.0, 1, false)
    @test_throws ArgumentError measuremap(model, NaN, state, 1, false)
    @test_throws ArgumentError measuremap(model, 1.0, state, 7, false)
    @test_throws ArgumentError measure_basismap(model, 1.0, ntuple(_ -> UInt8(0),6), 1, false)
end
