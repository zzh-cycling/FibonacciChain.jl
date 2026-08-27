using FibonacciChain
using Test

@testset "basis.jl" begin
    include("./test_Basis.jl")
end

@testset "Observable.jl" begin
    include("./test_Observable.jl")
end

@testset "AnyonLadder.jl" begin
    include("./test_AnyonLadder.jl")
end

@testset "Measurement.jl" begin
    include("./test_Measurement.jl")
end

@testset "MPSMeasurement.jl" begin
    include("./test_MPSMeasurement.jl")
end

@testset "HybridEvolution.jl" begin
    include("./test_HybridEvolution.jl")
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

@testset "OBF" begin
    include("./test_OBF.jl")
end

@testset "Heisenberg" begin
    include("./test_Heisenberg.jl")
end

@testset "test_topo_symmetry.jl" begin
    include("./test_topo_symmetry.jl")
end
