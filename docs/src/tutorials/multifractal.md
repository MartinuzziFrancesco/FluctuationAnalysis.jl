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

A [`MFDFAResult`](@ref) exposes the moment-resolved quantities through stable
accessors:

```@example mf
generalized_hurst(result)      # h(q)
```

```@example mf
mass_exponents(result)         # τ(q) = q·h(q) − 1
```

```@example mf
singularity_strengths(result)  # α
```

```@example mf
singularity_spectrum(result)   # f(α)
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
maximum(generalized_hurst(result)) - minimum(generalized_hurst(result))
```

A multiplicative binomial cascade, by contrast, is strongly multifractal. Here
the spread of ``h(q)`` is wide:

```@example mf
function binomial_cascade(levels, multiplier, rng)
    measure = [1.0]
    for _ in 1:levels
        refined = Vector{Float64}(undef, 2 * length(measure))
        for (index, value) in pairs(measure)
            left, right = if rand(rng) < 0.5
                (multiplier, 1 - multiplier)
            else
                (1 - multiplier, multiplier)
            end
            refined[2index - 1] = value * left
            refined[2index]     = value * right
        end
        measure = refined
    end
    return measure
end

cascade = binomial_cascade(14, 0.3, MersenneTwister(7))
cascade_result = mfdfa(cascade; q_values=collect(-4.0:0.5:4.0))
maximum(generalized_hurst(cascade_result)) - minimum(generalized_hurst(cascade_result))
```

The width of the singularity spectrum ``f(\alpha)`` measures the strength of the
multifractality:

```@example mf
strengths = singularity_strengths(cascade_result)
maximum(strengths) - minimum(strengths)
```

## Multifractal DMA

[`mfdma`](@ref) takes the same `q_values` and returns an [`MFDMAResult`](@ref)
with the same accessors, using a moving-average detrending whose position is set
by `theta`:

```@example mf
mfdma(cascade; q_values=collect(-4.0:0.5:4.0), theta=0.0) |> scaling_exponent
```
```
