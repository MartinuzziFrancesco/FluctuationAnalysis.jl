@doc doc"""
    MFDMAResult

Result of a multifractal detrending moving average analysis, returned by
[`mfdma`](@ref).

# Public interface

Use [`analysis_scales`](@ref), [`fluctuation_values`](@ref),
[`analysis_method`](@ref), [`fit_results`](@ref), [`moment_orders`](@ref),
[`generalized_hurst`](@ref), [`mass_exponents`](@ref),
[`singularity_strengths`](@ref), [`singularity_spectrum`](@ref), and
[`scaling_exponent`](@ref) to inspect the result.

The accessors preserve `float(eltype(series))`: moment orders are cast to the
data type, so a `Float32` series yields `Float32` values. Concrete fields are
implementation details and are not part of the stable public interface.
"""
@concrete struct MFDMAResult <: AbstractFluctuationResult
    q_values
    scales::Vector{Int}
    fluctuations
    moving_average::MovingAverage
    fits
    generalized_hurst
    mass_exponents
    singularity_strengths
    singularity_spectrum
end

Base.show(stream::IO, result::MFDMAResult) = __show_multifractal_summary(stream, result)

@doc doc"""
    scaling_exponent(result::Union{MFDFAResult, MFDMAResult}) -> Real

Generalized Hurst exponent at ``q = 2``, the standard DFA/DMA scaling exponent.

# Arguments

- `result::Union{MFDFAResult, MFDMAResult}`: a multifractal result whose moment
  orders include ``q = 2``.

# Returns

- the exponent ``h(2)``, in the result's value type.

# Throws

- `ArgumentError`: if ``q = 2`` is not among the moment orders.
"""
function scaling_exponent(result::Union{MFDFAResult, MFDMAResult})
    index = findfirst(Base.Fix2(isapprox, 2), result.q_values)
    index === nothing &&
        throw(ArgumentError("scaling_exponent requires q = 2 among the q values"))
    return result.generalized_hurst[index]
end

@doc doc"""
    mfdma(series; kwargs...) -> MFDMAResult

Multifractal detrending moving average (MFDMA) analysis of a one-dimensional time
`series`, following Gu & Zhou (2010).

A global moving-average trend is removed from the integrated profile at each scale,
the residual is split into disjoint segments, and the segment mean squares are
combined into q-order fluctuations for every order in `q_values`. The log-log slope
of each gives a generalized Hurst exponent ``h(q)``, from which the mass exponents
``\tau(q) = q\,h(q) - 1`` and, by a Legendre transform, the singularity strengths
``\alpha`` and spectrum ``f(\alpha)`` are derived. The exact ``q=0`` case uses the
logarithmic limit. At ``q = 2`` the result matches [`dma`](@ref) exactly.

# Arguments

- `series::AbstractVector{<:Real}`: the time series; must have at least 8 points.

# Keywords

- `q_values::AbstractVector{<:Real} = collect(-5.0:0.5:5.0)`: moment orders;
  sorted and deduplicated, with at least two distinct values required.
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
  `(lower, upper)` scale bounds restricting each log-log fit.

# Returns

- [`MFDMAResult`](@ref): the q values, scales, fluctuation matrix, moving average,
  fits, and the derived ``h(q)``, ``\tau(q)``, ``\alpha``, and ``f(\alpha)``.

# Throws

- `ArgumentError`: if `series` has fewer than 8 points or is constant, if fewer
  than two distinct finite `q` values are given, if a scale yields a zero-variance
  segment while a nonpositive moment is requested, or if a fitted fluctuation is
  not positive and finite.
"""
function mfdma(
        series::AbstractVector{<:Real};
        q_values::AbstractVector{<:Real} = collect(-5.0:0.5:5.0),
        theta::Real = 0.0,
        moving_average::MovingAverage = MovingAverage(theta),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        demean::Bool = false,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for MFDMA"))
    __require_nonconstant_series(series, "MFDMA")
    sorted_q = sort(unique(float.(q_values)))
    all(isfinite, sorted_q) || throw(ArgumentError("q values must be finite"))
    length(sorted_q) >= 2 ||
        throw(ArgumentError("need at least two distinct q values for MFDMA"))

    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = __mfdma_fluctuations(profile, scales, sorted_q, moving_average)
    # Carry the fluctuation value type through every derived quantity.
    q_values_typed = eltype(fluctuations).(sorted_q)
    fits = __fit_generalized_hurst(
        scales, fluctuations, q_values_typed; fitrange = fitrange
    )
    hurst_values = [fit.exponent for fit in fits]
    mass_exponent_values = __compute_mass_exponents(q_values_typed, hurst_values)
    strengths, spectrum = __compute_singularity_spectrum(
        q_values_typed, mass_exponent_values
    )

    return MFDMAResult(
        q_values_typed,
        scales,
        fluctuations,
        moving_average,
        fits,
        hurst_values,
        mass_exponent_values,
        strengths,
        spectrum,
    )
end
