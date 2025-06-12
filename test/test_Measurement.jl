using FibonacciChain
using Test
using LinearAlgebra

@testset "measure_basismap" begin
    N = 3
    ϕ = (1 + √5) / 2
    T = BitStr{N, Int}
    state = T(0b000)
    idx = 2
    pbc = false
    τ =0.0
    cstτ = 1/2(exp(τ)+1)/√(exp(2τ)+1)
    sign = :p
    basis0 = [T(0b000), T(0b001), T(0b010), T(0b100), T(0b101)]

    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, cstτ)
    @test output[2] == (T(bit"001"), cstτ)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, cstτ)
    @test output[4] == (T(bit"100"), cstτ)
    @test output[5] == (T(bit"101"), cstτ)

    τ = 0.0
    sign = :m
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ, cstτ)
    @test output[2] == (T(bit"001"), cstτ)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ, cstτ)
    @test output[4] == (T(bit"100"), cstτ)
    @test output[5] == (T(bit"101"), cstτ)

    τ = 1.0
    sign = :p
    cstτ = 1/2(exp(τ)+1)/√(exp(2τ)+1)
    coef = (exp(τ)-1)/2√(exp(2τ)+1)
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), cstτ-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), cstτ-2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    τ = 1.0
    sign = :m
    cstτ = 1/2(exp(τ)+1)/√(exp(2τ)+1)
    coef = (1-exp(τ))/2√(exp(2τ)+1)
    output = measure_basismap.(T, τ, basis0, idx, sign, pbc)
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"010"), cstτ+coef*(1-2ϕ^(-1)), cstτ-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), cstτ-2*coef*ϕ^(-3/2))
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] == (T(bit"101"), cstτ-coef)

    idx = 1
    output = measure_basismap.(T, τ, basis0, idx, sign) # pbc is true by default
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"100"), cstτ+coef*(1-2ϕ^(-1)), cstτ-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), cstτ+coef)
    @test output[3] == (T(bit"010"), cstτ+coef)
    @test output[4] == (T(bit"100"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), cstτ-2*coef*ϕ^(-3/2))
    @test output[5] === nothing

    idx = 3
    output = measure_basismap.(T, τ, basis0, idx, sign) # pbc is true by default
    @test length(output) == length(basis0)
    @test output[1] == (T(bit"000"), T(bit"001"), cstτ+coef*(1-2ϕ^(-1)), cstτ-2*coef*ϕ^(-3/2))
    @test output[2] == (T(bit"001"), T(bit"000"), cstτ+coef*(2ϕ^(-1)-1), cstτ-2*coef*ϕ^(-3/2))
    @test output[3] == (T(bit"010"), cstτ+coef)
    @test output[4] == (T(bit"100"), cstτ+coef)
    @test output[5] === nothing
end

@testset "measure_matrix" begin
    N = 3
    T = BitStr{N, Int}
    ϕ = (1 + √5) / 2


    τ = 1.0
    idx = 2
    cstτ = (exp(1)+1)/2√(exp(2)+1)
    coef = (exp(1)-1)/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 cstτ-2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        cstτ-2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef]
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, :p, false)
    @test Mpobc == expected_matrix 

    coef = (1-exp(1))/2√(exp(2)+1)
    expected_matrix = [
        cstτ+coef*(1-2ϕ^(-1)) 0.0 cstτ-2*coef*ϕ^(-3/2) 0.0 0.0;
        0.0 cstτ+coef 0.0 0.0 0.0;
        cstτ-2*coef*ϕ^(-3/2)  0.0 cstτ+coef*(2ϕ^(-1)-1) 0.0 0.0;
        0.0 0.0 0.0 cstτ+coef 0.0;
        0.0 0.0 0.0 0.0 cstτ-coef
    ]
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, :m, false)
    @test Mmobc == expected_matrix
    @test Mpobc^2+Mmobc^2 == I(5) 

    expected_matrix = []
    Mppbc = FibonacciChain.measure_matrix(T, τ, idx, :p) # pbc
    @test Mppbc == expected_matrix 
    Mmpbc = FibonacciChain.measure_matrix(T, τ, idx, :m) # pbc
    @test Mmpbc == expected_matrix    
    @test Mppbc^2+Mmpbc^2 == I(T) "The sum of Mppbc^2 and Mmpbc^2 should be the identity matrix, but got $(Mppbc^2 + Mmpbc^2)"

    # Test with a different τ
    τ = 0.0          
    expected_matrix = []
    Mpobc = FibonacciChain.measure_matrix(T, τ, idx, :p, pbc)
    @test Mpobc == expected_matrix 
    expected_matrix = []
    Mmobc = FibonacciChain.measure_matrix(T, τ, idx, :m, pbc)
    @test Mmobc == expected_matrix 
    @test Mpobc^2+Mmobc^2 == I(T) "The sum of Mpobc^2 and Mmobc^2 should be the identity matrix, but got $(Mpobc^2 + Mmobc^2)"

    expected_matrix = []
    Mppbc = FibonacciChain.measure_matrix(T, τ, idx, :p) # pbc
    @test Mppbc == expected_matrix 
    expected_matrix = []
    Mmpbc = FibonacciChain.measure_matrix(T, τ, idx, :m) # pbc
    @test Mmpbc == expected_matrix 
  
    @test Mppbc^2+Mmpbc^2 == I(T) "The sum of Mppbc^2 and Mmpbc^2 should be the identity matrix, but got $(Mppbc^2 + Mmpbc^2)"
    # Test with a different idx'
    idx = 1
    sign = :p
    expected_matrix = []

end

@testset "measuremap" begin
    N = 3
    T = BitStr{N, Int}
    τ = 1.0
    idx = 2
    sign = :p
    pbc = true
    state = T(0b000)
    expected_output = T(0b000)
    output = measuremap(T, τ, idx, sign, state, pbc)        
    @test output == expected_output "Output should be $(expected_output), but got $(output)"
    @test output isa T "Output should be of type $(T), but got $(typeof(output))"
    @test output == T(0b000) "Output should be T(0b000), but got $(output)"
    τ = 0.0
    sign = :n
    expected_output = T(0b000)
    output = measuremap(T, τ, idx, sign, state, pbc)
    @test output == expected_output "Output should be $(expected_output), but got $(output)"
    @test output isa T "Output should be of type $(T), but got $(typeof(output))"
    @test output == T(0b000) "Output should be T(0b000), but got $(output)"
    # Test with a different state
    state = T(0b001)
    expected_output = T(0b001)
    output = measuremap(T, τ, idx, sign, state, pbc)
    @test output == expected_output "Output should be $(expected_output), but got $(output)"
    @test output isa T "Output should be of type $(T), but got $(typeof(output))"
    @test output == T(0b001) "Output should be T(0b001), but got $(output)"
end