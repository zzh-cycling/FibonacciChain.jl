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

@testset "MPSMeasurement.jl" begin
    include("./test_MPSMeasurement.jl")
end

@testset "FiboSparse.jl" begin
    include("./test_FiboSparse.jl")
end

@testset "Ising" begin
    include("./test_Ising.jl")
end

@testset "test_Reference.jl" begin
    include("./test_Reference.jl")
end