using FibonacciChain
using Test

@testset "basis.jl" begin
    include("./test_Basis.jl")
end

@testset "basis.jl" begin
    include("./test_Observable.jl")
end

@testset "ladderFibo.jl" begin
    include("./test_LadderFibo.jl")
end