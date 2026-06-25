"""
    dfa(series; kwargs...) -> DFAResult

Detrended fluctuation analysis of a one-dimensional time `series`.

The series is converted into its integrated profile, the profile is split into
segments at each scale, a local trend is removed by the detrender, and the
root-mean-square fluctuation is regressed against scale in log-log coordinates.
The slope is the DFA scaling exponent, read with [`scaling_exponent`](@ref).

# Arguments

- `series::AbstractVector{<:Real}`: the time series; must have at least 8 points.

# Keywords

- `order::Integer = 1`: polynomial detrending order used to build the default
  detrender.
- `detrender::AbstractDetrender = PolynomialDetrender(order)`: detrender removing
  the local trend; overrides `order` when given.
- `scales::AbstractVector{<:Integer} = logarithmic_scales(length(series))`:
  window sizes to evaluate.
- `demean::Bool = true`: subtract the mean before integrating the profile.
- `overlap::Bool = false`: use overlapping sliding segments instead of disjoint
  ones.
- `bidirectional::Bool = true`: also segment from the end of the profile when not
  overlapping, so trailing samples are covered.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`:
  `(lower, upper)` scale bounds restricting the log-log fit.

# Returns

- [`DFAResult`](@ref): the scales, fluctuation values, detrender, and log-log fit.

# Throws

- `ArgumentError`: if `series` has fewer than 8 points, or if a scale is too small
  for the chosen detrender.
"""
function dfa(
        series::AbstractVector{<:Real};
        order::Integer = 1,
        detrender::AbstractDetrender = PolynomialDetrender(order),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        demean::Bool = true,
        overlap::Bool = false,
        bidirectional::Bool = true,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for DFA"))

    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = fluctuation_curve(
        profile, scales, detrender; overlap = overlap, bidirectional = bidirectional
    )
    fit = loglog_fit(scales, fluctuations; fitrange = fitrange)

    return DFAResult(scales, fluctuations, detrender, fit)
end
