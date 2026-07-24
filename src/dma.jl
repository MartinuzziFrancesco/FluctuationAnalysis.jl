"""
    DMAResult

Result of a detrending moving average analysis, returned by [`dma`](@ref).

# Public interface

Use [`analysis_scales`](@ref), [`fluctuation_values`](@ref),
[`analysis_method`](@ref), [`fit_results`](@ref), and
[`scaling_exponent`](@ref) to inspect the result. The fluctuation values and fit
preserve `float(eltype(series))`.

Concrete fields are implementation details and are not part of the stable public
interface.
"""
@concrete struct DMAResult <: AbstractFluctuationResult
    scales::Vector{Int}
    fluctuations
    moving_average::MovingAverage
    fit
end

Base.show(stream::IO, result::DMAResult) = __show_scaling_summary(stream, result)

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
- `demean::Bool = false`: subtract the mean before integrating the profile. Gu &
  Zhou (2010) eq. (1) uses the raw cumulative sum, so this defaults to `false`;
  demeaning leaves an offset under one-sided (`theta = 0` or `1`) moving averages
  and biases their exponents, so it should stay `false` for those.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`:
  `(lower, upper)` scale bounds restricting the log-log fit.

# Returns

- [`DMAResult`](@ref): the scales, fluctuation values, moving average, and fit.

# Throws

- `ArgumentError`: if `series` has fewer than 8 points or is constant, if a scale
  cannot form a residual segment, or if the fitted fluctuations are not positive
  and finite.
"""
function dma(
        series::AbstractVector{<:Real};
        theta::Real = 0.0,
        moving_average::MovingAverage = MovingAverage(theta),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        demean::Bool = false,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for DMA"))
    __require_nonconstant_series(series, "DMA")

    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = __dma_fluctuation_curve(profile, scales, moving_average)
    fit = loglog_fit(scales, fluctuations; fitrange = fitrange)

    return DMAResult(scales, fluctuations, moving_average, fit)
end
