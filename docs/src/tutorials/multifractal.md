```@meta
CurrentModule = FluctuationAnalysis
```

# Multifractal analysis

A monofractal signal is described by a single scaling exponent. Many real signals
are *multifractal*: their scaling depends on the magnitude of the fluctuations
being weighted. [`mfdfa`](@ref) resolves this by raising the segment variances to
a range of moment orders ``q`` before fitting, yielding a generalized Hurst
exponent ``h(q)`` for each order. [`mfdma`](@ref) does the same with a
moving-average detrending (see [Detrending moving average](@ref)).

## Running MFDFA

```@example mf
using FluctuationAnalysis
using Random

series = randn(MersenneTwister(1), 30_000)
result = mfdfa(series; q_values=collect(-4.0:1.0:4.0))
result
```

A [`MFDFAResult`](@ref) exposes the moment-resolved quantities directly:

```@example mf
result.generalized_hurst      # h(q)
```

```@example mf
result.mass_exponents         # τ(q) = q·h(q) − 1
```

```@example mf
result.singularity_strengths  # α
```

```@example mf
result.singularity_spectrum   # f(α)
```

At ``q = 2`` the analysis coincides with [`dfa`](@ref), so
[`scaling_exponent`](@ref) returns ``h(2)``:

```@example mf
scaling_exponent(result)
```

## Monofractal versus multifractal

For white noise the generalized Hurst exponents are nearly flat — a narrow spread
of ``h(q)`` is the signature of a monofractal signal:

```@example mf
maximum(result.generalized_hurst) - minimum(result.generalized_hurst)
```

A multiplicative binomial cascade, by contrast, is strongly multifractal. Here
the spread of ``h(q)`` is wide:

```@example mf
function binomial_cascade(levels, multiplier, rng)
    measure = [1.0]
    for _ in 1:levels
        refined = Vector{Float64}(undef, 2 * length(measure))
        for (index, value) in pairs(measure)
            left, right = rand(rng) < 0.5 ? (multiplier, 1 - multiplier) : (1 - multiplier, multiplier)
            refined[2index - 1] = value * left
            refined[2index]     = value * right
        end
        measure = refined
    end
    return measure
end

cascade = binomial_cascade(14, 0.3, MersenneTwister(7))
cascade_result = mfdfa(cascade; q_values=collect(-4.0:0.5:4.0))
maximum(cascade_result.generalized_hurst) - minimum(cascade_result.generalized_hurst)
```

The width of the singularity spectrum ``f(\alpha)`` measures the strength of the
multifractality:

```@example mf
maximum(cascade_result.singularity_strengths) - minimum(cascade_result.singularity_strengths)
```

## Multifractal DMA

[`mfdma`](@ref) takes the same `q_values` and returns an [`MFDMAResult`](@ref)
with the same fields, using a moving-average detrending whose position is set by
`theta`:

```@example mf
mfdma(cascade; q_values=collect(-4.0:0.5:4.0), theta=0.0) |> scaling_exponent
```
```
