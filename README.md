# FluctuationAnalysis [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl/dev/) [![Build Status](https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MartinuzziFrancesco/FluctuationAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/MartinuzziFrancesco/FluctuationAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MartinuzziFrancesco/FluctuationAnalysis.jl) [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) [![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for fluctuation and scaling analysis of time series. It provides
detrended fluctuation analysis (`dfa`) and its multifractal (`mfdfa`) and
cross-correlation (`dcca`) extensions, the detrending moving average family
(`dma`, `mfdma`), and Hurst exponent estimation (`hurst`). The package is built
from reusable components (profile construction, scale generation, segmentation,
detrending, fluctuation computation, and log-log fitting) that every method
shares.

See the [documentation](https://MartinuzziFrancesco.github.io/FluctuationAnalysis.jl/dev/)
for tutorials and the full API reference.

## Usage

```julia
using FluctuationAnalysis

series = randn(10_000)
result = dfa(series)

result.fit.exponent   # scaling exponent
result.scales         # window sizes
result.fluctuations   # fluctuation function values
```

The detrending order, the scales, and the segmentation can all be configured:

```julia
result = dfa(series; order=2, scales=logarithmic_scales(length(series)),
            overlap=false, bidirectional=true)
```

A custom detrender can be supplied through the `AbstractDetrender` interface:

```julia
result = dfa(series; detrender=PolynomialDetrender(3))
```

## Multifractal DFA

`mfdfa` generalizes DFA over a range of moment orders `q`, returning the
generalized Hurst exponents `h(q)`, the mass exponents `τ(q)`, and the
singularity spectrum `f(α)`:

```julia
result = mfdfa(series; q_values=collect(-5.0:0.5:5.0))

result.generalized_hurst      # h(q)
result.mass_exponents         # τ(q)
result.singularity_strengths  # α
result.singularity_spectrum   # f(α)
scaling_exponent(result)      # h(2), the standard DFA exponent
```

## Detrending moving average

`dma` is an alternative to `dfa` that removes the local trend with a moving
average instead of a polynomial fit; `mfdma` is its multifractal version. The
window position is set by `theta` (`0` backward, `0.5` centered, `1` forward):

```julia
result = dma(series; theta=0.0)
scaling_exponent(result)         # DMA scaling exponent

mfdma(series; q_values=collect(-5.0:0.5:5.0))   # multifractal version
```

## Detrended cross-correlation

`dcca` analyzes two equal-length series with the same segmentation and detrending
as DFA, returning the cross-correlation exponent, the detrended cross-covariance,
and the DCCA cross-correlation coefficient `ρ_DCCA(s)`:

```julia
result = dcca(first_series, second_series)

scaling_exponent(result)   # cross-correlation exponent λ
result.correlation         # ρ_DCCA(s), in [-1, 1]
result.covariances         # detrended cross-covariance F²_DCCA(s)
```

When the two series are identical the analysis reduces to `dfa`.

## Hurst exponent

The Hurst exponent can be estimated either from the DFA scaling exponent (the
default) or from rescaled range (R/S) analysis:

```julia
hurst(series)                                # DFA-based estimate
hurst(series, RescaledRangeHurst())          # rescaled range estimate

result = hurst(series)
hurst_exponent(result)   # the estimated exponent
```
