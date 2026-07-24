function __q_order_fluctuation(variances::AbstractVector{<:Real}, order_q::Real)
    if iszero(order_q)
        return exp(mean(log.(variances)) / 2)
    end
    order = convert(float(eltype(variances)), order_q)
    order == 2 && return sqrt(mean(variances))
    any(iszero, variances) && return mean(variances .^ (order / 2))^(1 / order)
    log_variances = log.(variances)
    log_center = mean(log_variances)
    centered_moment = mean(exp.((order / 2) .* (log_variances .- log_center)))
    return exp(log_center / 2 + log(centered_moment) / order)
end

function __q_order_fluctuations!(
        row::AbstractVector{<:Real},
        variances::AbstractVector{<:Real},
        q_values::AbstractVector{<:Real},
        scale::Integer,
    )
    all(>=(0), variances) ||
        throw(ArgumentError("scale $scale produced a negative segment variance"))
    if any(<=(0), q_values) && any(iszero, variances)
        throw(
            ArgumentError(
                "scale $scale produced a zero-variance segment; nonpositive " *
                    "multifractal moments are undefined",
            ),
        )
    end
    for (q_index, order_q) in pairs(q_values)
        row[q_index] = __q_order_fluctuation(variances, order_q)
    end
    return row
end

function __mfdfa_fluctuations(
        profile::AbstractVector{T},
        scales::AbstractVector{<:Integer},
        q_values::AbstractVector{<:Real},
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    ) where {T <: Real}
    smallest_allowed = __minimum_segment_length(detrender)
    fluctuations = zeros(T, length(scales), length(q_values))
    for (scale_index, scale) in pairs(scales)
        scale >= smallest_allowed ||
            throw(ArgumentError("scale $scale is too small for the chosen detrender"))
        variances = __segment_variances(
            profile, scale, detrender; overlap = overlap, bidirectional = bidirectional
        )
        __q_order_fluctuations!(
            view(fluctuations, scale_index, :), variances, q_values, scale
        )
    end
    return fluctuations
end

function __fit_generalized_hurst(
        scales::AbstractVector{<:Integer},
        fluctuations::AbstractMatrix{<:Real},
        q_values::AbstractVector{<:Real};
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    return [
        loglog_fit(scales, view(fluctuations, :, q_index); fitrange = fitrange) for
            q_index in eachindex(q_values)
    ]
end

function __compute_mass_exponents(
        q_values::AbstractVector{<:Real}, hurst_values::AbstractVector{<:Real}
    )
    return q_values .* hurst_values .- 1
end

function __central_difference(nodes::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    node_count = length(nodes)
    length(values) == node_count ||
        throw(ArgumentError("nodes and values must have equal length"))
    node_count >= 2 || throw(ArgumentError("need at least two points to differentiate"))
    for index in 2:node_count
        nodes[index] > nodes[index - 1] ||
            throw(ArgumentError("nodes must be strictly increasing"))
    end
    value_type = promote_type(float(eltype(nodes)), float(eltype(values)))
    derivative = zeros(value_type, node_count)
    derivative[1] = (values[2] - values[1]) / (nodes[2] - nodes[1])
    derivative[node_count] =
        (values[node_count] - values[node_count - 1]) /
        (nodes[node_count] - nodes[node_count - 1])
    for index in 2:(node_count - 1)
        previous_node = nodes[index - 1]
        current_node = nodes[index]
        next_node = nodes[index + 1]
        previous_weight =
            (current_node - next_node) /
            ((previous_node - current_node) * (previous_node - next_node))
        current_weight =
            (2current_node - previous_node - next_node) /
            ((current_node - previous_node) * (current_node - next_node))
        next_weight =
            (current_node - previous_node) /
            ((next_node - previous_node) * (next_node - current_node))
        derivative[index] =
            previous_weight * values[index - 1] + current_weight * values[index] +
            next_weight * values[index + 1]
    end
    return derivative
end

function __compute_singularity_spectrum(
        q_values::AbstractVector{<:Real}, mass_exponent_values::AbstractVector{<:Real}
    )
    strengths = __central_difference(q_values, mass_exponent_values)
    spectrum = q_values .* strengths .- mass_exponent_values
    return strengths, spectrum
end

@doc doc"""
    MFDFAResult

Result of a multifractal detrended fluctuation analysis, returned by
[`mfdfa`](@ref).

# Public interface

Use [`analysis_scales`](@ref), [`fluctuation_values`](@ref),
[`analysis_method`](@ref), [`fit_results`](@ref), [`moment_orders`](@ref),
[`generalized_hurst`](@ref), [`mass_exponents`](@ref),
[`singularity_strengths`](@ref), [`singularity_spectrum`](@ref), and
[`scaling_exponent`](@ref) to inspect the result.

The accessors preserve `float(eltype(series))`: moment orders are cast to the
data type, so a `Float32` series yields `Float32` values regardless of the
`q_values` type. Concrete fields are implementation details and are not part of
the stable public interface.
"""
@concrete struct MFDFAResult <: AbstractFluctuationResult
    q_values
    scales::Vector{Int}
    fluctuations
    detrender
    fits
    generalized_hurst
    mass_exponents
    singularity_strengths
    singularity_spectrum
end

Base.show(stream::IO, result::MFDFAResult) = __show_multifractal_summary(stream, result)

@doc doc"""
    mfdfa(series; kwargs...) -> MFDFAResult

Multifractal detrended fluctuation analysis of a one-dimensional time `series`.

The detrended segment variances are combined into q-order fluctuations for every
order in `q_values`; the log-log slope of each gives a generalized Hurst exponent
``h(q)``, from which the mass exponents ``\tau(q) = q\,h(q) - 1`` and, by a
Legendre transform, the singularity strengths ``\alpha`` and spectrum
``f(\alpha)`` are derived. The exact ``q=0`` case uses the logarithmic limit. At
``q = 2`` the analysis reduces exactly to [`dfa`](@ref).

# Arguments

- `series::AbstractVector{<:Real}`: the time series; must have at least 8 points.

# Keywords

- `q_values::AbstractVector{<:Real} = collect(-5.0:0.5:5.0)`: moment orders;
  sorted and deduplicated, with at least two distinct values required.
- `order::Integer = 1`: polynomial detrending order for the default detrender.
- `detrender::AbstractDetrender = PolynomialDetrender(order)`: detrender removing
  the local trend; overrides `order` when given.
- `scales::AbstractVector{<:Integer} = logarithmic_scales(length(series))`:
  window sizes to evaluate.
- `demean::Bool = true`: subtract the mean before integrating the profile.
- `overlap::Bool = false`: use overlapping sliding segments instead of disjoint
  ones.
- `bidirectional::Bool = true`: also segment from the end of the profile when not
  overlapping.
- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`:
  `(lower, upper)` scale bounds restricting each log-log fit.

# Returns

- [`MFDFAResult`](@ref): the q values, scales, fluctuation matrix, fits, and the
  derived ``h(q)``, ``\tau(q)``, ``\alpha``, and ``f(\alpha)``.

# Throws

- `ArgumentError`: if `series` has fewer than 8 points or is constant, if fewer
  than two distinct finite `q` values are given, if a scale yields a zero-variance
  segment while a nonpositive moment is requested, or if a fitted fluctuation is
  not positive and finite.
"""
function mfdfa(
        series::AbstractVector{<:Real};
        q_values::AbstractVector{<:Real} = collect(-5.0:0.5:5.0),
        order::Integer = 1,
        detrender::AbstractDetrender = PolynomialDetrender(order),
        scales::AbstractVector{<:Integer} = logarithmic_scales(length(series)),
        demean::Bool = true,
        overlap::Bool = false,
        bidirectional::Bool = true,
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(series) >= 8 || throw(ArgumentError("series is too short for MFDFA"))
    __require_nonconstant_series(series, "MFDFA")
    sorted_q = sort(unique(float.(q_values)))
    all(isfinite, sorted_q) || throw(ArgumentError("q values must be finite"))
    length(sorted_q) >= 2 ||
        throw(ArgumentError("need at least two distinct q values for MFDFA"))
    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = __mfdfa_fluctuations(
        profile,
        scales,
        sorted_q,
        detrender;
        overlap = overlap,
        bidirectional = bidirectional,
    )
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

    return MFDFAResult(
        q_values_typed,
        scales,
        fluctuations,
        detrender,
        fits,
        hurst_values,
        mass_exponent_values,
        strengths,
        spectrum,
    )
end
