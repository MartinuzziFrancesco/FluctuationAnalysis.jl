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

- [`dfa`](@ref): detrended fluctuation analysis, following
  [Peng et al. (1994)](https://doi.org/10.1103/PhysRevE.49.1685), with
  bidirectional segmentation from
  [Kantelhardt et al. (2002)](https://doi.org/10.1016/S0378-4371(02)01383-3).
- [`mfdfa`](@ref): multifractal DFA, following
  [Kantelhardt et al. (2002)](https://doi.org/10.1016/S0378-4371(02)01383-3).
- [`dcca`](@ref): detrended cross-correlation analysis, following
  [Podobnik and Stanley (2008)](https://doi.org/10.1103/PhysRevLett.100.084102),
  with the coefficient from
  [Zebende (2011)](https://doi.org/10.1016/j.physa.2010.10.022).
- [`dma`](@ref): detrending moving average, the ``q=2`` MFDMA case, following
  [Gu and Zhou (2010)](https://doi.org/10.1103/PhysRevE.82.011136).
- [`mfdma`](@ref): multifractal DMA, following
  [Gu and Zhou (2010)](https://doi.org/10.1103/PhysRevE.82.011136).
- [`hurst`](@ref): DFA-based or rescaled-range Hurst estimation, following
  [Hurst (1951)](https://doi.org/10.1061/TACEAT.0006518) and
  [Peng et al. (1994)](https://doi.org/10.1103/PhysRevE.49.1685).

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
