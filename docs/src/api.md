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

## Detrender interface

```@docs
AbstractDetrender
PolynomialDetrender
detrend
minimum_segment_length
```

## Hurst estimator interface

```@docs
AbstractHurstEstimator
DetrendedFluctuationHurst
RescaledRangeHurst
hurst_statistic_curve
```
```
