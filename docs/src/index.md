```@meta
CurrentModule = FluctuationAnalysis
```

# FluctuationAnalysis.jl

[FluctuationAnalysis.jl](https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl)
provides fluctuation and scaling analysis of one-dimensional time series. It
implements detrended fluctuation analysis and its multifractal,
cross-correlation, and moving-average relatives, together with Hurst exponent
estimation, on top of a small set of reusable building blocks.

## Methods

| Function          | Method                                              | Primary reference                                                              |
|:----------------- |:--------------------------------------------------- |:----------------------------------------------------------------------------- |
| [`dfa`](@ref)     | Detrended fluctuation analysis (DFA)                | [Peng et al. (1994)](https://doi.org/10.1103/PhysRevE.49.1685)                 |
| [`mfdfa`](@ref)   | Multifractal DFA                                    | [Kantelhardt et al. (2002)](https://doi.org/10.1016/S0378-4371(02)01383-3)     |
| [`dcca`](@ref)    | Detrended cross-correlation analysis (DCCA)         | [Podobnik & Stanley (2008)](https://doi.org/10.1103/PhysRevLett.100.084102)    |
| [`dma`](@ref)     | Detrending moving average (DMA)                     | [Alessio et al. (2002)](https://doi.org/10.1140/epjb/e20020150)                |
| [`mfdma`](@ref)   | Multifractal DMA                                    | [Gu & Zhou (2010)](https://doi.org/10.1103/PhysRevE.82.011136)                 |
| [`hurst`](@ref)   | Hurst exponent (DFA-based or rescaled range)        | [Hurst (1951)](https://doi.org/10.1061/TACEAT.0006518); [Peng et al. (1994)](https://doi.org/10.1103/PhysRevE.49.1685) |

Every method shares the same pipeline — [`integrated_profile`](@ref),
[`logarithmic_scales`](@ref), segmentation, detrending through the
[`AbstractDetrender`](@ref) interface, and a [`loglog_fit`](@ref) — and returns a
rich result object whose scaling exponent is read with [`scaling_exponent`](@ref).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl")
```

## Quick start

```@example index
using FluctuationAnalysis
using Random

series = randn(MersenneTwister(1), 10_000)
result = dfa(series)
scaling_exponent(result)
```

## Where to go next

- New to the package? Start with [Getting started](@ref), which walks through a
  full DFA workflow and the options shared by every method.
- For the moment-resolved multifractal spectrum, see
  [Multifractal analysis](@ref).
- The [API reference](@ref) lists every exported function and type.
```
