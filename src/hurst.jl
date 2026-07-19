"""
    AbstractHurstEstimator

Supertype for Hurst exponent estimators.

Each estimator defines how a scale-dependent statistic is built from a series;
the Hurst exponent is then the slope of that statistic against scale in log-log
coordinates.

See also [`DetrendedFluctuationHurst`](@ref) and [`RescaledRangeHurst`](@ref).
"""
abstract type AbstractHurstEstimator end

"""
    DetrendedFluctuationHurst(detrender) -> DetrendedFluctuationHurst
    DetrendedFluctuationHurst(; order = 1) -> DetrendedFluctuationHurst

Estimate the Hurst exponent from the detrended fluctuation analysis scaling
exponent. For a stationary signal the DFA exponent equals the Hurst exponent.
This is the default estimator used by [`hurst`](@ref).

# Arguments

- `detrender::AbstractDetrender`: detrender supplying the local trend model.

# Keywords

- `order::Integer = 1`: polynomial degree used to build the default detrender,
  for the keyword constructor.
"""
@concrete struct DetrendedFluctuationHurst <: AbstractHurstEstimator
    detrender
end

function DetrendedFluctuationHurst(; order::Integer = 1)
    return DetrendedFluctuationHurst(PolynomialDetrender(order))
end

"""
    RescaledRangeHurst() -> RescaledRangeHurst

Estimate the Hurst exponent from rescaled range (R/S) analysis. The rescaled
range averaged over non-overlapping segments grows with segment length as a power
law whose exponent is the Hurst exponent.
"""
@concrete struct RescaledRangeHurst <: AbstractHurstEstimator end

"""
    HurstResult

Result of a Hurst exponent estimation, returned by [`hurst`](@ref).

# Fields

- `estimator::AbstractHurstEstimator`: the estimator used.
- `scales::Vector{Int}`: segment lengths at which the statistic was evaluated.
- `statistic::Vector{T}`: the scale-dependent statistic, one per scale.
- `fit::LogLogFit{T}`: the fit whose `exponent` is the Hurst exponent.

The statistic and fit share the value type `float(eltype(series))`, preserving the
input precision.
"""
@concrete struct HurstResult <: AbstractFluctuationResult
    estimator
    scales::Vector{Int}
    statistic
    fit
end

Base.show(stream::IO, result::HurstResult) = __show_scaling_summary(stream, result)

function __rescaled_range(segment::AbstractVector{<:Real})
    segment_mean = mean(segment)
    cumulative_deviation = cumsum(segment .- segment_mean)
    range_width = maximum(cumulative_deviation) - minimum(cumulative_deviation)
    spread = std(segment; corrected = false, mean = segment_mean)
    spread == 0 && return zero(range_width / oneunit(spread))
    return range_width / spread
end

function __mean_rescaled_range(series::AbstractVector{<:Real}, scale::Integer)
    segments = __segment_views(series, scale; overlap = false, bidirectional = false)
    return mean(__rescaled_range, segments)
end

function __hurst_statistic_curve(
        estimator::DetrendedFluctuationHurst,
        series::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
    )
    profile = integrated_profile(series; demean = true)
    return __fluctuation_curve(profile, scales, estimator.detrender)
end

function __hurst_statistic_curve(
        ::RescaledRangeHurst, series::AbstractVector{<:Real}, scales::AbstractVector{<:Integer}
    )
    values = float.(collect(series))
    statistic = similar(values, length(scales))
    for (index, scale) in pairs(scales)
        statistic[index] = __mean_rescaled_range(values, scale)
    end
    return statistic
end

"""
    hurst(series, estimator = DetrendedFluctuationHurst(); kwargs...) -> HurstResult

Estimate the Hurst exponent of a one-dimensional time `series`.

The chosen `estimator` builds a scale-dependent statistic which is regressed
against scale in log-log coordinates; the slope is the Hurst exponent.

# Arguments

- `series::AbstractVector{<:Real}`: the time series; must have at least 8 points.
- `estimator::AbstractHurstEstimator = DetrendedFluctuationHurst()`: use
  [`DetrendedFluctuationHurst`](@ref) for the DFA-based estimate or
  [`RescaledRangeHurst`](@ref) for rescaled range (R/S) analysis.

# Keywords

- `scales::AbstractVector{<:Integer} = logarithmic_scales(length(series))`:
  segment lengths to evaluate.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`: restrict the
  fit to a `(lower, upper)` band of scales.

# Returns

- [`HurstResult`](@ref): the estimator, scales, statistic, and the log-log fit
  whose exponent is the Hurst exponent (also available via
  [`hurst_exponent`](@ref)).

# Throws

- `ArgumentError`: if `series` has fewer than 8 points.
"""
function hurst(
        series::AbstractVector{<:Real},
        estimator::AbstractHurstEstimator = DetrendedFluctuationHurst();
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for Hurst estimation"))
    scales = Int.(collect(scales))
    statistic = __hurst_statistic_curve(estimator, series, scales)
    fit = loglog_fit(scales, statistic; fitrange = fitrange)
    return HurstResult(estimator, scales, statistic, fit)
end

"""
    hurst_exponent(result::HurstResult) -> AbstractFloat

Estimated Hurst exponent held by a [`HurstResult`](@ref); equivalent to
[`scaling_exponent`](@ref) on the same result.

# Arguments

- `result::HurstResult`: a result returned by [`hurst`](@ref).

# Returns

- the estimated Hurst exponent, in the result's value type.
"""
hurst_exponent(result::HurstResult) = result.fit.exponent
