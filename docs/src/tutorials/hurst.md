```@meta
CurrentModule = FluctuationAnalysis
```

# Hurst exponent

The Hurst exponent ``H`` summarizes the long-range dependence of a series:
``H = 0.5`` for an uncorrelated process, ``H > 0.5`` for persistent (trending)
behaviour, and ``H < 0.5`` for anti-persistent behaviour. [`hurst`](@ref)
estimates it through a choice of [`AbstractHurstEstimator`](@ref).

## The default estimator

By default [`hurst`](@ref) uses [`DetrendedFluctuationHurst`](@ref), which reads
the exponent off a DFA fit — for a stationary signal the DFA exponent equals the
Hurst exponent:

```@example hurst
using FluctuationAnalysis
using Random

series = randn(MersenneTwister(1), 10_000)
result = hurst(series)
hurst_exponent(result)
```

[`hurst`](@ref) returns a [`HurstResult`](@ref); the exponent is available both
as [`hurst_exponent`](@ref) and through the generic [`scaling_exponent`](@ref).

```@example hurst
scaling_exponent(result) == hurst_exponent(result)
```

## Rescaled range analysis

Passing [`RescaledRangeHurst`](@ref) switches to classic R/S analysis, where the
rescaled range averaged over segments grows as a power law in the segment length:

```@example hurst
hurst_exponent(hurst(series, RescaledRangeHurst()))
```

R/S has a known positive small-sample bias, so on finite white noise it usually
sits slightly above ``0.5``.

## Choosing the detrending order

The DFA-based estimator accepts a detrending order, mirroring [`dfa`](@ref):

```@example hurst
hurst_exponent(hurst(series, DetrendedFluctuationHurst(; order=2)))
```

## Persistent signals

A cumulative sum is strongly persistent, and the estimate rises accordingly:

```@example hurst
walk = cumsum(randn(MersenneTwister(2), 10_000))
hurst_exponent(hurst(walk))
```
```
