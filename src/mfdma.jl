@doc doc"""
    MFDMAResult

Result of a multifractal detrending moving average analysis, returned by
[`mfdma`](@ref).

# Fields

- `q_values::Vector{T}`: moment orders, sorted ascending.
- `scales::Vector{Int}`: window sizes at which the fluctuations were evaluated.
- `fluctuations::Matrix{T}`: q-order fluctuations, indexed `[scale_index, q_index]`.
- `moving_average::MovingAverage`: the moving-average specification used to detrend.
- `fits::Vector{LogLogFit{T}}`: the log-log fit of each `q` value.
- `generalized_hurst::Vector{T}`: generalized Hurst exponents ``h(q)``.
- `mass_exponents::Vector{T}`: mass scaling exponents ``\tau(q)``.
- `singularity_strengths::Vector{T}`: singularity strengths ``\alpha``.
- `singularity_spectrum::Vector{T}`: singularity spectrum ``f(\alpha)``.

The result's value type is `float(eltype(series))`: the moment orders are cast to
the data type, so a `Float32` series yields a `Float32` result.
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

function Base.show(stream::IO, result::MFDMAResult)
    print(
        stream,
        "MFDMAResult(q=",
        length(result.q_values),
        " in [",
        round(minimum(result.q_values); digits = 2),
        ", ",
        round(maximum(result.q_values); digits = 2),
        "], scales=",
        length(result.scales),
        ")",
    )
    return nothing
end

@doc doc"""
    scaling_exponent(result::MFDMAResult) -> AbstractFloat

Generalized Hurst exponent at ``q = 2``, the standard DMA scaling exponent.

# Arguments

- `result::MFDMAResult`: a multifractal result whose moment orders include
  ``q = 2``.

# Returns

- the exponent ``h(2)``, in the result's value type.

# Throws

- `ArgumentError`: if ``q = 2`` is not among the moment orders.
"""
function scaling_exponent(result::MFDMAResult)
    index = findfirst(order_q -> isapprox(order_q, 2), result.q_values)
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
``\alpha`` and spectrum ``f(\alpha)`` are derived. At ``q = 2`` the result matches
[`dma`](@ref).

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
- `demean::Bool = true`: subtract the mean before integrating the profile.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`:
  `(lower, upper)` scale bounds restricting each log-log fit.

# Returns

- [`MFDMAResult`](@ref): the q values, scales, fluctuation matrix, moving average,
  fits, and the derived ``h(q)``, ``\tau(q)``, ``\alpha``, and ``f(\alpha)``.

# Throws

- `ArgumentError`: if `series` has fewer than 8 points, if fewer than two distinct
  `q` values are given, or if a scale yields a zero-variance segment.
"""
function mfdma(
        series::AbstractVector{<:Real};
        q_values::AbstractVector{<:Real} = collect(-5.0:0.5:5.0),
        theta::Real = 0.0,
        moving_average::MovingAverage = MovingAverage(theta),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        demean::Bool = true,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for MFDMA"))
    sorted_q = sort(unique(float.(q_values)))
    length(sorted_q) >= 2 ||
        throw(ArgumentError("need at least two distinct q values for MFDMA"))

    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = mfdma_fluctuations(profile, scales, sorted_q, moving_average)
    # Carry the fluctuation value type through every derived quantity.
    q_values_typed = eltype(fluctuations).(sorted_q)
    fits = fit_generalized_hurst(scales, fluctuations, q_values_typed; fitrange = fitrange)
    hurst_values = [fit.exponent for fit in fits]
    mass_exponent_values = compute_mass_exponents(q_values_typed, hurst_values)
    strengths, spectrum = compute_singularity_spectrum(q_values_typed, mass_exponent_values)

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
