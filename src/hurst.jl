"""
    AbstractHurstEstimator

Supertype for Hurst exponent estimators.

Each estimator defines how a scale-dependent statistic is built from a series;
the Hurst exponent is then the slope of that statistic against scale in log-log
coordinates.

# Interface

A concrete subtype `E <: AbstractHurstEstimator` must implement
`hurst_statistic_curve(estimator::E, series, scales)`, returning the statistic at
each scale.

See also [`DetrendedFluctuationHurst`](@ref), [`RescaledRangeHurst`](@ref), and
[`hurst_statistic_curve`](@ref).
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
struct RescaledRangeHurst <: AbstractHurstEstimator end

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

function Base.show(stream::IO, result::HurstResult)
    print(
        stream,
        "HurstResult(exponent=",
        round(result.fit.exponent; digits = 4),
        ", scales=",
        length(result.scales),
        ", rsquared=",
        round(result.fit.rsquared; digits = 4),
        ")",
    )
    return nothing
end

"""
    rescaled_range(segment)

Rescaled range (R/S) of a `segment`: the range of its cumulative mean-removed
deviations divided by its population standard deviation. A segment with zero
spread has a rescaled range of zero.
"""
function rescaled_range(segment::AbstractVector{<:Real})
    segment_mean = mean(segment)
    cumulative_deviation = cumsum(segment .- segment_mean)
    range_width = maximum(cumulative_deviation) - minimum(cumulative_deviation)
    spread = std(segment; corrected = false, mean = segment_mean)
    spread == 0 && return zero(range_width / oneunit(spread))
    return range_width / spread
end

"""
    mean_rescaled_range(series, scale)

Mean rescaled range over the non-overlapping segments of `series` of length
`scale`.
"""
function mean_rescaled_range(series::AbstractVector{<:Real}, scale::Integer)
    segments = segment_views(series, scale; overlap = false, bidirectional = false)
    return mean(rescaled_range, segments)
end

"""
    hurst_statistic_curve(estimator::AbstractHurstEstimator, series, scales) -> Vector{<:AbstractFloat}

Scale-dependent statistic whose log-log slope estimates the Hurst exponent.

Part of the [`AbstractHurstEstimator`](@ref) interface, with methods for
[`DetrendedFluctuationHurst`](@ref) (the DFA fluctuation curve) and
[`RescaledRangeHurst`](@ref) (the mean rescaled range).

# Arguments

- `estimator::AbstractHurstEstimator`: the estimator selecting the statistic.
- `series::AbstractVector{<:Real}`: the time series.
- `scales::AbstractVector{<:Integer}`: segment lengths to evaluate.

# Returns

- `Vector{<:AbstractFloat}`: the statistic at each scale.
"""
function hurst_statistic_curve(
        estimator::DetrendedFluctuationHurst,
        series::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
    )
    profile = integrated_profile(series; demean = true)
    return fluctuation_curve(profile, scales, estimator.detrender)
end

function hurst_statistic_curve(
        ::RescaledRangeHurst, series::AbstractVector{<:Real}, scales::AbstractVector{<:Integer}
    )
    values = float.(collect(series))
    statistic = similar(values, length(scales))
    for (index, scale) in pairs(scales)
        statistic[index] = mean_rescaled_range(values, scale)
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
    statistic = hurst_statistic_curve(estimator, series, scales)
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
