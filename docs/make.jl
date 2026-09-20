using FibonacciChain
using Documenter

DocMeta.setdocmeta!(
    FibonacciChain,
    :DocTestSetup,
    :(using FibonacciChain, Random, LinearAlgebra, BitBasis, ITensorMPS, ITensors);
    recursive = true,
)

makedocs(;
    modules = [FibonacciChain],
    doctest = true,
    authors = "Zhaohui Zhi",
    sitename = "FibonacciChain.jl",
    format = Documenter.HTML(;
        canonical = "https://zzh-cycling.github.io/FibonacciChain.jl",
        assets = String[],
        edit_link = nothing,
        # The index intentionally contains the complete generated API reference.
        size_threshold = 500 * 2^10,
        size_threshold_warn = 300 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Basis Functions" => "basis.md",
            "Ising Anyon Chain" => "ising_anyons.md",
            "Observables" => "observables.md",
            "Measurements" => "measurements.md",
            "Hybrid Evolution" => "hybrid_evolution.md",
            "MPS Methods" => "mps.md",
            "Examples" => "examples.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/zzh-cycling/FibonacciChain.jl",
    devbranch = "main",
    forcepush = true,
)
