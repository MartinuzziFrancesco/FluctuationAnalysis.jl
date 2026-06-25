using FluctuationAnalysis
using Documenter

DocMeta.setdocmeta!(
    FluctuationAnalysis, :DocTestSetup, :(using FluctuationAnalysis); recursive = true
)

makedocs(;
    modules = [FluctuationAnalysis],
    authors = "Francesco Martinuzzi",
    sitename = "FluctuationAnalysis.jl",
    format = Documenter.HTML(;
        canonical = "https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl",
        edit_link = "main",
        assets = String[],
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting started" => "tutorials/getting_started.md",
            "Detrending moving average" => "tutorials/moving_average.md",
            "Multifractal analysis" => "tutorials/multifractal.md",
            "Cross-correlation" => "tutorials/cross_correlation.md",
            "Hurst exponent" => "tutorials/hurst.md",
        ],
        "API reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/MartinuzziFrancesco/FluctuationAnalysis.jl", devbranch = "main")
