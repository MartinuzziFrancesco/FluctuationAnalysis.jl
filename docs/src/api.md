```@meta
CurrentModule = FluctuationAnalysis
```

# API reference

```@docs
FluctuationAnalysis.FluctuationAnalysis
```

```@index
```

## Detrended fluctuation analysis

```@docs
dfa
DFAResult
```

## Multifractal DFA

```@docs
mfdfa
MFDFAResult
```

## Detrended cross-correlation analysis

```@docs
dcca
DCCAResult
```

## Moving-average methods

```@docs
dma
DMAResult
mfdma
MFDMAResult
MovingAverage
```

## Hurst estimation

```@docs
hurst
hurst_exponent
HurstResult
```

## Shared building blocks

```@docs
integrated_profile
logarithmic_scales
loglog_fit
LogLogFit
scaling_exponent
AbstractFluctuationResult
```

## Result accessors

Use these accessors instead of depending on the concrete field layout of result
types.

```@docs
analysis_scales
fluctuation_values
fit_results
analysis_method
moment_orders
generalized_hurst
mass_exponents
singularity_strengths
singularity_spectrum
dcca_covariances
dcca_correlation
dcca_marginal_fluctuations
hurst_statistic
```

## Detrender interface

```@docs
AbstractDetrender
PolynomialDetrender
detrend
```

## Hurst estimator interface

```@docs
AbstractHurstEstimator
DetrendedFluctuationHurst
RescaledRangeHurst
```
