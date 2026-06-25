```@meta
CurrentModule = FluctuationAnalysis
```

# Cross-correlation

Detrended cross-correlation analysis (DCCA) extends DFA to a *pair* of series. It
measures how their detrended fluctuations move together as a function of scale,
returning both a cross-correlation exponent and a scale-resolved correlation
coefficient. The entry point is [`dcca`](@ref).

## Two coupled series

We build two series that share a common random component plus independent noise,
so they are genuinely correlated:

```@example dcca
using FluctuationAnalysis
using Random

rng = MersenneTwister(1)
common = cumsum(randn(rng, 10_000))
first_series  = common .+ randn(rng, 10_000)
second_series = common .+ randn(rng, 10_000)

result = dcca(first_series, second_series)
result
```

## The correlation coefficient

The DCCA cross-correlation coefficient ``\rho_{DCCA}(s)`` lies in ``[-1, 1]`` and
is defined at every scale, even where the cross-correlation exponent is not:

```@example dcca
result.correlation
```

Because the two series share a strong common component, the coefficient is high
across scales:

```@example dcca
using Statistics
mean(result.correlation)
```

For independent series it collapses toward zero:

```@example dcca
independent = dcca(randn(MersenneTwister(2), 10_000), randn(MersenneTwister(3), 10_000))
mean(independent.correlation)
```

## The cross-correlation exponent

[`scaling_exponent`](@ref) returns the slope of the detrended cross-covariance
over the scales where it is positive:

```@example dcca
scaling_exponent(result)
```

A [`DCCAResult`](@ref) also stores the signed covariance and the single-series
fluctuations used to form the coefficient:

```@example dcca
(cov=result.covariances[1], fx=result.first_fluctuations[1], fy=result.second_fluctuations[1])
```

## Reduction to DFA

When both inputs are the same series, DCCA reduces exactly to [`dfa`](@ref): the
cross-covariance becomes the variance and ``\rho_{DCCA}(s) \equiv 1``.

```@example dcca
self = dcca(first_series, first_series)
(exponent=scaling_exponent(self), rho=mean(self.correlation))
```
```
