```@meta
CurrentModule = FluctuationAnalysis
```

# Detrending moving average

The detrending moving average (DMA) method, introduced by Gu & Zhou (2010),
replaces the per-segment polynomial fit of [`dfa`](@ref) with a *global* moving
average of the integrated profile. [`dma`](@ref) is its single-exponent form and
[`mfdma`](@ref) its multifractal generalization.

## Running DMA

```@example dma
using FluctuationAnalysis
using Random

series = randn(MersenneTwister(1), 10_000)
result = dma(series)
scaling_exponent(result)
```

[`dma`](@ref) returns a [`DMAResult`](@ref) with the same shape as a
[`DFAResult`](@ref):

```@example dma
result.fluctuations
```

## The window position

The moving-average window can sit behind, around, or ahead of each point. The
position is the `theta` keyword, captured by the [`MovingAverage`](@ref) type:

- `theta = 0` — *backward* (past points only), the default and Gu & Zhou's most
  accurate variant for multifractal measures;
- `theta = 0.5` — *centered*, the only variant that removes a linear trend
  exactly;
- `theta = 1` — *forward* (future points only).

```@example dma
[theta => scaling_exponent(dma(series; theta=theta)) for theta in (0.0, 0.5, 1.0)]
```

For all three, white noise scales with an exponent near ``0.5``.

## DMA versus DFA

DMA estimates a trend with a simple moving average rather than a least-squares
polynomial, so its detrending is weaker than higher-order DFA and it tends to
underestimate large exponents. The two agree closely on uncorrelated signals:

```@example dma
(dfa=scaling_exponent(dfa(series)), dma=scaling_exponent(dma(series)))
```

## Multifractal DMA

[`mfdma`](@ref) combines the moving-average residuals across moment orders `q`,
exactly as [`mfdfa`](@ref) does for the polynomial residuals. See
[Multifractal analysis](@ref) for the shape of the output; at `q = 2` it reduces
to [`dma`](@ref):

```@example dma
m = mfdma(series; q_values=[2.0, 4.0])
isapprox(scaling_exponent(m), scaling_exponent(dma(series)); atol=1e-10)
```
```
