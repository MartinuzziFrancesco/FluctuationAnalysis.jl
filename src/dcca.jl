"""
    cross_covariance_curve(first_profile, second_profile, scales, detrender;
                           overlap=false, bidirectional=true)

Detrended cross-covariance ``F^2_{DCCA}(s)`` evaluated over every entry of
`scales`, averaged over the aligned segment pairs of the two profiles. The
values are signed.
"""
function cross_covariance_curve(
        first_profile::AbstractVector{<:Real},
        second_profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    smallest_allowed = minimum_segment_length(detrender)
    value_type = promote_type(float(eltype(first_profile)), float(eltype(second_profile)))
    covariances = Vector{value_type}(undef, length(scales))
    for (index, scale) in pairs(scales)
        scale >= smallest_allowed ||
            throw(ArgumentError("scale $scale is too small for the chosen detrender"))
        covariances[index] = mean(
            segment_covariances(
                first_profile,
                second_profile,
                scale,
                detrender;
                overlap = overlap,
                bidirectional = bidirectional,
            ),
        )
    end
    return covariances
end

@doc doc"""
    DCCAResult

Result of a detrended cross-correlation analysis of two time series, returned by
[`dcca`](@ref).

# Fields

- `scales::Vector{Int}`: window sizes at which the analysis was evaluated.
- `covariances::Vector{T}`: detrended cross-covariance ``F^2_{DCCA}(s)`` (signed).
- `cross_fluctuations::Vector{T}`: signed cross fluctuation
  ``\mathrm{sign}(F^2)\sqrt{|F^2|}``.
- `first_fluctuations::Vector{T}`: DFA fluctuation of the first series.
- `second_fluctuations::Vector{T}`: DFA fluctuation of the second series.
- `correlation::Vector{T}`: DCCA cross-correlation coefficient ``\rho_{DCCA}(s)``
  in ``[-1, 1]``.
- `detrender::AbstractDetrender`: detrender used to remove local trends.
- `fit::LogLogFit{T}`: the fit whose `exponent` is the cross-correlation exponent,
  taken over the scales with positive detrended covariance.

The stored vectors and fit share the value type, the promotion of the two series'
floated element types.
"""
@concrete struct DCCAResult <: AbstractFluctuationResult
    scales::Vector{Int}
    covariances
    cross_fluctuations
    first_fluctuations
    second_fluctuations
    correlation
    detrender
    fit
end

function Base.show(stream::IO, result::DCCAResult)
    print(
        stream,
        "DCCAResult(exponent=",
        round(result.fit.exponent; digits = 4),
        ", scales=",
        length(result.scales),
        ", mean_rho=",
        round(mean(result.correlation); digits = 4),
        ")",
    )
    return nothing
end

"""
    fit_cross_scaling(scales, covariances, cross_fluctuations; fitrange=nothing)

Fit the cross-correlation exponent from the scales whose detrended covariance is
positive and which lie within `fitrange`.
"""
function fit_cross_scaling(
        scales::AbstractVector{<:Integer},
        covariances::AbstractVector{<:Real},
        cross_fluctuations::AbstractVector{<:Real};
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    selection = scale_selection(scales, fitrange) .& (covariances .> 0)
    count(selection) >= 2 || throw(
        ArgumentError(
            "fewer than two scales with positive detrended covariance in the fit range; " *
                "a cross-correlation exponent is not defined, inspect the `correlation` field instead",
        ),
    )
    return loglog_fit(scales[selection], cross_fluctuations[selection])
end

@doc doc"""
    dcca(first_series, second_series; kwargs...) -> DCCAResult

Detrended cross-correlation analysis (DCCA) of two equal-length time series.

Both series are converted into integrated profiles and split into segments at each
scale; the detrended covariance of each aligned segment pair is averaged into
``F^2_{DCCA}(s)``. The cross-correlation exponent is the log-log slope of
``\sqrt{F^2_{DCCA}(s)}`` over the scales where the covariance is positive. The
DCCA cross-correlation coefficient ``\rho_{DCCA}(s) \in [-1, 1]`` is returned for
every scale. When the two series are identical the analysis reduces to
[`dfa`](@ref).

# Arguments

- `first_series::AbstractVector{<:Real}`: the first time series.
- `second_series::AbstractVector{<:Real}`: the second time series; must be the
  same length as `first_series` (at least 8 points).

# Keywords

- `order::Integer = 1`: polynomial detrending order for the default detrender.
- `detrender::AbstractDetrender = PolynomialDetrender(order)`: detrender removing
  the local trend; overrides `order` when given.
- `scales::AbstractVector{<:Integer} = logarithmic_scales(length(first_series))`:
  window sizes to evaluate.
- `demean::Bool = true`: subtract the mean before integrating each profile.
- `overlap::Bool = false`: use overlapping sliding segments instead of disjoint
  ones.
- `bidirectional::Bool = true`: also segment from the end of the profile when not
  overlapping.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`:
  `(lower, upper)` scale bounds restricting the log-log fit.

# Returns

- [`DCCAResult`](@ref): the scales, signed covariance, cross fluctuations, both
  single-series fluctuations, the ``\rho_{DCCA}(s)`` coefficient, and the fit.

# Throws

- `ArgumentError`: if the series differ in length or have fewer than 8 points, or
  if fewer than two scales have positive detrended covariance in the fit range
  (in which case inspect the `correlation` field instead of the exponent).
"""
function dcca(
        first_series::AbstractVector{<:Real},
        second_series::AbstractVector{<:Real};
        order::Integer = 1,
        detrender::AbstractDetrender = PolynomialDetrender(order),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(first_series)),
        demean::Bool = true,
        overlap::Bool = false,
        bidirectional::Bool = true,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(first_series) == length(second_series) ||
        throw(ArgumentError("the two series must have equal length"))
    length(first_series) >= 8 || throw(ArgumentError("series is too short for DCCA"))

    scales = Int.(collect(scales))
    first_profile = integrated_profile(first_series; demean = demean)
    second_profile = integrated_profile(second_series; demean = demean)

    covariances = cross_covariance_curve(
        first_profile,
        second_profile,
        scales,
        detrender;
        overlap = overlap,
        bidirectional = bidirectional,
    )
    # Unify both single-series fluctuation curves to the covariance value type so
    # every stored vector shares one element type (a no-op when types already match).
    value_type = eltype(covariances)
    first_fluctuations = convert(
        Vector{value_type},
        fluctuation_curve(
            first_profile, scales, detrender; overlap = overlap, bidirectional = bidirectional
        ),
    )
    second_fluctuations = convert(
        Vector{value_type},
        fluctuation_curve(
            second_profile, scales, detrender; overlap = overlap, bidirectional = bidirectional
        ),
    )

    cross_fluctuations = sign.(covariances) .* sqrt.(abs.(covariances))
    correlation = covariances ./ (first_fluctuations .* second_fluctuations)
    fit = fit_cross_scaling(scales, covariances, cross_fluctuations; fitrange = fitrange)

    return DCCAResult(
        scales,
        covariances,
        cross_fluctuations,
        first_fluctuations,
        second_fluctuations,
        correlation,
        detrender,
        fit,
    )
end
