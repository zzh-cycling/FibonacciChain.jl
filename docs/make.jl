using FibonacciChain
using Documenter

DocMeta.setdocmeta!(FibonacciChain, :DocTestSetup, :(using FibonacciChain, Random, LinearAlgebra, BitBasis, ITensorMPS, ITensors); recursive=true)

makedocs(;
    modules=[FibonacciChain],
    doctest=true,
    authors="Zhaohui Zhi",
    sitename="FibonacciChain.jl",
    format=Documenter.HTML(;
        canonical="https://zzh-cycling.github.io/FibonacciChain.jl",
        assets=String[],
        edit_link=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Basis Functions" => "basis.md",
            "Observables" => "observables.md", 
            "Measurements" => "measurements.md",
            "MPS Methods" => "mps.md",
            "Examples" => "examples.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/zzh-cycling/FibonacciChain.jl",
    devbranch="main",
    forcepush=true,
)
