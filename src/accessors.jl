const __MULTIFRACTAL_RESULTS = Union{MFDFAResult, MFDMAResult}

"""
    analysis_scales(result::AbstractFluctuationResult) -> Vector{Int}

Scales at which an analysis evaluated its scale-dependent quantity.

# Arguments

- `result::AbstractFluctuationResult`: a result returned by an exported analysis.

# Returns

- the ordered integer scales stored by `result`.
"""
analysis_scales(result::AbstractFluctuationResult) = result.scales

@doc doc"""
    fluctuation_values(result) -> AbstractArray

Fluctuation values evaluated at each scale.

For DFA and DMA results this is a vector indexed by scale. For multifractal
results it is a matrix indexed by `(scale_index, q_index)`. For DCCA results it
is the signed cross-fluctuation ``\operatorname{sign}(F^2)\sqrt{|F^2|}``.

# Arguments

- `result::Union{DFAResult,MFDFAResult,DCCAResult,DMAResult,MFDMAResult}`: a
  fluctuation-analysis result.

# Returns

- the fluctuation vector or matrix stored by `result`.
"""
fluctuation_values(result::DFAResult) = result.fluctuations
fluctuation_values(result::MFDFAResult) = result.fluctuations
fluctuation_values(result::DCCAResult) = result.cross_fluctuations
fluctuation_values(result::DMAResult) = result.fluctuations
fluctuation_values(result::MFDMAResult) = result.fluctuations

"""
    fit_results(result::AbstractFluctuationResult) -> Tuple{Vararg{LogLogFit}}

Log-log fits used to estimate the result's scaling exponents.

Scalar analyses return a one-element tuple. Multifractal analyses return one fit
per moment order, in the same order as [`moment_orders`](@ref).

# Arguments

- `result::AbstractFluctuationResult`: a result returned by an exported analysis.

# Returns

- a tuple of [`LogLogFit`](@ref) objects.
"""
fit_results(result::DFAResult) = (result.fit,)
fit_results(result::DCCAResult) = (result.fit,)
fit_results(result::DMAResult) = (result.fit,)
fit_results(result::HurstResult) = (result.fit,)
fit_results(result::__MULTIFRACTAL_RESULTS) = Tuple(result.fits)

"""
    analysis_method(result::AbstractFluctuationResult) -> Any

Method specification used to construct an analysis result.

This is the detrender for DFA, MFDFA, and DCCA; the moving-average specification
for DMA and MFDMA; and the Hurst estimator for a Hurst result.

# Arguments

- `result::AbstractFluctuationResult`: a result returned by an exported analysis.

# Returns

- the result's `AbstractDetrender`, `MovingAverage`, or `AbstractHurstEstimator`.
"""
analysis_method(result::Union{DFAResult, MFDFAResult, DCCAResult}) = result.detrender
analysis_method(result::Union{DMAResult, MFDMAResult}) = result.moving_average
analysis_method(result::HurstResult) = result.estimator

"""
    moment_orders(result) -> AbstractVector

Moment orders evaluated by a multifractal analysis.

# Arguments

- `result::$(__MULTIFRACTAL_RESULTS)`: a multifractal result.

# Returns

- the sorted, distinct moment orders stored by `result`.
"""
moment_orders(result::__MULTIFRACTAL_RESULTS) = result.q_values

@doc doc"""
    generalized_hurst(result) -> AbstractVector

Generalized Hurst exponents ``h(q)`` from a multifractal analysis.

# Arguments

- `result::$(__MULTIFRACTAL_RESULTS)`: a multifractal result.

# Returns

- one generalized Hurst exponent per value in [`moment_orders`](@ref).
"""
generalized_hurst(result::__MULTIFRACTAL_RESULTS) = result.generalized_hurst

@doc doc"""
    mass_exponents(result) -> AbstractVector

Mass scaling exponents ``\tau(q)`` from a multifractal analysis.

# Arguments

- `result::$(__MULTIFRACTAL_RESULTS)`: a multifractal result.

# Returns

- one mass exponent per value in [`moment_orders`](@ref).
"""
mass_exponents(result::__MULTIFRACTAL_RESULTS) = result.mass_exponents

@doc doc"""
    singularity_strengths(result) -> AbstractVector

Singularity strengths ``\alpha`` from a multifractal analysis.

# Arguments

- `result::$(__MULTIFRACTAL_RESULTS)`: a multifractal result.

# Returns

- the singularity-strength vector derived from the mass exponents.
"""
singularity_strengths(result::__MULTIFRACTAL_RESULTS) = result.singularity_strengths

@doc doc"""
    singularity_spectrum(result) -> AbstractVector

Singularity spectrum ``f(\alpha)`` from a multifractal analysis.

# Arguments

- `result::$(__MULTIFRACTAL_RESULTS)`: a multifractal result.

# Returns

- the singularity-spectrum vector paired with [`singularity_strengths`](@ref).
"""
singularity_spectrum(result::__MULTIFRACTAL_RESULTS) = result.singularity_spectrum

@doc doc"""
    dcca_covariances(result::DCCAResult) -> AbstractVector

Signed detrended cross-covariances ``F^2_{DCCA}(s)``.

# Arguments

- `result::DCCAResult`: a result returned by [`dcca`](@ref).

# Returns

- one signed detrended cross-covariance per scale.
"""
dcca_covariances(result::DCCAResult) = result.covariances

@doc doc"""
    dcca_correlation(result::DCCAResult) -> AbstractVector

DCCA cross-correlation coefficients ``\rho_{DCCA}(s)``.

# Arguments

- `result::DCCAResult`: a result returned by [`dcca`](@ref).

# Returns

- one coefficient in `[-1, 1]` per scale.
"""
dcca_correlation(result::DCCAResult) = result.correlation

"""
    dcca_marginal_fluctuations(result::DCCAResult) -> NamedTuple

DFA fluctuation curves of the two DCCA input series.

# Arguments

- `result::DCCAResult`: a result returned by [`dcca`](@ref).

# Returns

- a named tuple `(first, second)` containing the two marginal fluctuation vectors.
"""
function dcca_marginal_fluctuations(result::DCCAResult)
    return (first = result.first_fluctuations, second = result.second_fluctuations)
end

"""
    hurst_statistic(result::HurstResult) -> AbstractVector

Scale-dependent statistic fitted by a Hurst analysis.

# Arguments

- `result::HurstResult`: a result returned by [`hurst`](@ref).

# Returns

- one DFA fluctuation or rescaled-range statistic per scale, according to the
  result's estimator.
"""
hurst_statistic(result::HurstResult) = result.statistic
