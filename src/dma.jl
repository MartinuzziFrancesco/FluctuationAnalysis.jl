"""
    DMAResult

Result of a detrending moving average analysis, returned by [`dma`](@ref).

# Fields

- `scales::Vector{Int}`: window sizes at which the fluctuation was evaluated.
- `fluctuations::Vector{T}`: fluctuation function values, one per scale.
- `moving_average::MovingAverage`: the moving-average specification used to detrend.
- `fit::LogLogFit{T}`: the fit holding the scaling exponent and metadata.

The fluctuations and fit share the value type `float(eltype(series))`, preserving
the input precision.
"""
@concrete struct DMAResult <: AbstractFluctuationResult
    scales::Vector{Int}
    fluctuations
    moving_average::MovingAverage
    fit
end

function Base.show(stream::IO, result::DMAResult)
    print(
        stream,
        "DMAResult(exponent=",
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
    dma(series; kwargs...) -> DMAResult

Detrending moving average (DMA) analysis of a one-dimensional time `series`,
following Gu & Zhou (2010).

The series is converted into its integrated profile, a global moving-average trend
is removed at each scale, the residual is split into disjoint segments, and the
root-mean-square fluctuation is regressed against scale in log-log coordinates.
Its slope is the DMA scaling exponent, equal to [`mfdma`](@ref) at `q = 2`.

# Arguments

- `series::AbstractVector{<:Real}`: the time series; must have at least 8 points.

# Keywords

- `theta::Real = 0.0`: moving-average position in `[0, 1]` (`0` backward, `0.5`
  centered, `1` forward); used to build the default moving average.
- `moving_average::MovingAverage = MovingAverage(theta)`: the moving-average
  specification; overrides `theta` when given.
- `scales::AbstractVector{<:Integer} = logarithmic_scales(length(series))`:
  window sizes to evaluate.
- `demean::Bool = true`: subtract the mean before integrating the profile.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`:
  `(lower, upper)` scale bounds restricting the log-log fit.

# Returns

- [`DMAResult`](@ref): the scales, fluctuation values, moving average, and fit.

# Throws

- `ArgumentError`: if `series` has fewer than 8 points, or if a scale is too large
  to form a residual segment.
"""
function dma(
        series::AbstractVector{<:Real};
        theta::Real = 0.0,
        moving_average::MovingAverage = MovingAverage(theta),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        demean::Bool = true,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for DMA"))

    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = dma_fluctuation_curve(profile, scales, moving_average)
    fit = loglog_fit(scales, fluctuations; fitrange = fitrange)

    return DMAResult(scales, fluctuations, moving_average, fit)
end
