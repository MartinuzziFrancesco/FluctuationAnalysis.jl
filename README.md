# FluctuationAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl/dev/)
[![Build Status](https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MartinuzziFrancesco/FluctuationAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MartinuzziFrancesco/FluctuationAnalysis.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

FluctuationAnalysis.jl is a Julia package for fluctuation and scaling analysis of
time series. It is built from a small set of reusable, composable components —
profile construction, scale generation, segmentation, detrending, fluctuation
computation, and log-log fitting — that every method shares, so the analyses stay
consistent and new methods plug in cleanly. Each method returns a rich result
object carrying the scales, fluctuation values, fitted exponents, and diagnostics.

## Features

FluctuationAnalysis.jl provides functions to estimate scaling exponents and
characterize long-range correlations, multifractality, and cross-correlations in
time series. More specifically the software offers:

  - **Detrended fluctuation analysis** (`dfa`) with configurable polynomial
    detrending order, scales, overlapping/bidirectional segmentation, and a
    pluggable `AbstractDetrender` interface.
  - **Multifractal DFA** (`mfdfa`), returning the generalized Hurst exponents
    `h(q)`, the mass exponents `τ(q)`, the singularity strengths `α`, and the
    singularity spectrum `f(α)`.
  - **Detrended cross-correlation analysis** (`dcca`) for two series, returning
    the cross-correlation exponent `λ`, the detrended cross-covariance, and the
    DCCA cross-correlation coefficient `ρ_DCCA(s)`.
  - **Detrending moving average** family (`dma`, `mfdma`), an alternative to the
    polynomial detrending of DFA with a configurable window position `theta`.
  - **Hurst exponent estimation** (`hurst`), from either the DFA scaling exponent
    or rescaled range (R/S) analysis.
  - **Reusable building blocks** (`logarithmic_scales`, `PolynomialDetrender`,
    `MovingAverage`, fluctuation curves, `loglog_fit`) shared across every method
    and exposed for building custom analyses.

All methods are type-generic (`Float32`, `Float64`, `BigFloat`, …) and accept a
wide range of array inputs. See the
[documentation](https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl/dev/)
for tutorials and the full API reference.

## Installation

FluctuationAnalysis.jl is not registered yet. You can install it directly from
the repository using either of

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl")
```

or, from the Pkg REPL (press `]`):

```julia
pkg> add https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl
```

## Quick Example

The example below runs detrended fluctuation analysis on a white-noise series,
whose scaling exponent is close to `0.5`, and then extends the same data to its
multifractal spectrum:

```julia
using FluctuationAnalysis

series = randn(10_000)

result = dfa(series)
result.scales              # the window sizes
result.fluctuations        # the fluctuation function values
scaling_exponent(result)   # the DFA scaling exponent, ≈ 0.5

# tune the detrending order, scales, and segmentation
result = dfa(series; order = 2, overlap = false, bidirectional = true)

# multifractal spectrum over a range of moment orders q
mf = mfdfa(series; q_values = collect(-5.0:0.5:5.0))
mf.generalized_hurst       # h(q)
mf.singularity_strengths   # α
mf.singularity_spectrum    # f(α)
```
