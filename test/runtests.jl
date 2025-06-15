using FibonacciChain
using Test

@testset "basis.jl" begin
    include("./test_Basis.jl")
end

@testset "Observable.jl" begin
    include("./test_Observable.jl")
end

@testset "ladderFibo.jl" begin
    include("./test_LadderFibo.jl")
end

@testset "Measurement.jl" begin
    include("./test_Measurement.jl")
end