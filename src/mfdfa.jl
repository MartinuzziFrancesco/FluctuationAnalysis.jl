"""
    q_order_fluctuation(variances, order_q)

Combine the detrended segment `variances` at a scale into the q-order
fluctuation. For `order_q == 0` the geometric (logarithmic) average is used, as
prescribed for multifractal detrended fluctuation analysis. The moment order is
cast to the variance element type so the result keeps that type.
"""
function q_order_fluctuation(variances::AbstractVector{<:Real}, order_q::Real)
    if abs(float(order_q)) <= 16eps(Float64)
        return exp(mean(log.(variances)) / 2)
    end
    order = convert(float(eltype(variances)), order_q)
    return mean(variances .^ (order / 2))^(1 / order)
end

"""
    mfdfa_fluctuations(profile, scales, q_values, detrender; overlap=false, bidirectional=true)

Matrix of q-order fluctuations whose entry `[scale_index, q_index]` is the
fluctuation at `scales[scale_index]` for `q_values[q_index]`.
"""
function mfdfa_fluctuations(
        profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        q_values::AbstractVector{<:Real},
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    smallest_allowed = minimum_segment_length(detrender)
    value_type = float(eltype(profile))
    fluctuations = Matrix{value_type}(undef, length(scales), length(q_values))
    for (scale_index, scale) in pairs(scales)
        scale >= smallest_allowed ||
            throw(ArgumentError("scale $scale is too small for the chosen detrender"))
        variances = segment_variances(
            profile, scale, detrender; overlap = overlap, bidirectional = bidirectional
        )
        all(variance -> variance > 0, variances) || throw(
            ArgumentError(
                "scale $scale produced a zero-variance segment; multifractal moments are undefined",
            ),
        )
        for (q_index, order_q) in pairs(q_values)
            fluctuations[scale_index, q_index] = q_order_fluctuation(variances, order_q)
        end
    end
    return fluctuations
end

"""
    fit_generalized_hurst(scales, fluctuations, q_values; fitrange=nothing)

Fit each column of `fluctuations` against `scales` in log-log coordinates,
returning one [`LogLogFit`](@ref) per q value whose `exponent` is the
generalized Hurst exponent.
"""
function fit_generalized_hurst(
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

"""
    compute_mass_exponents(q_values, hurst_values)

Mass (Rényi) scaling exponents ``τ(q) = q·h(q) - 1`` for a one-dimensional
support.
"""
function compute_mass_exponents(
        q_values::AbstractVector{<:Real}, hurst_values::AbstractVector{<:Real}
    )
    return q_values .* hurst_values .- 1
end

"""
    central_difference(nodes, values)

Numerical derivative of `values` with respect to `nodes` using central
differences on the interior and one-sided differences at the endpoints.
"""
function central_difference(nodes::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    node_count = length(nodes)
    node_count >= 2 || throw(ArgumentError("need at least two points to differentiate"))
    derivative = Vector{promote_type(float(eltype(nodes)), float(eltype(values)))}(
        undef, node_count
    )
    derivative[1] = (values[2] - values[1]) / (nodes[2] - nodes[1])
    derivative[node_count] =
        (values[node_count] - values[node_count - 1]) /
        (nodes[node_count] - nodes[node_count - 1])
    for index in 2:(node_count - 1)
        derivative[index] =
            (values[index + 1] - values[index - 1]) / (nodes[index + 1] - nodes[index - 1])
    end
    return derivative
end

"""
    compute_singularity_spectrum(q_values, mass_exponent_values)

Singularity strengths ``α = dτ/dq`` and singularity spectrum
``f(α) = q·α - τ(q)`` obtained from the mass exponents by a Legendre transform.
"""
function compute_singularity_spectrum(
        q_values::AbstractVector{<:Real}, mass_exponent_values::AbstractVector{<:Real}
    )
    strengths = central_difference(q_values, mass_exponent_values)
    spectrum = q_values .* strengths .- mass_exponent_values
    return strengths, spectrum
end

@doc doc"""
    MFDFAResult

Result of a multifractal detrended fluctuation analysis, returned by
[`mfdfa`](@ref).

# Fields

- `q_values::Vector{T}`: moment orders, sorted ascending.
- `scales::Vector{Int}`: window sizes at which the fluctuations were evaluated.
- `fluctuations::Matrix{T}`: q-order fluctuations, indexed `[scale_index, q_index]`.
- `detrender::AbstractDetrender`: detrender used to remove local trends.
- `fits::Vector{LogLogFit{T}}`: the log-log fit of each `q` value.
- `generalized_hurst::Vector{T}`: generalized Hurst exponents ``h(q)``.
- `mass_exponents::Vector{T}`: mass scaling exponents ``\tau(q)``.
- `singularity_strengths::Vector{T}`: singularity strengths ``\alpha``.
- `singularity_spectrum::Vector{T}`: singularity spectrum ``f(\alpha)``.

The result's value type is `float(eltype(series))`: the moment orders are cast to
the data type, so a `Float32` series yields a `Float32` result regardless of the
`q_values` type.
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

function Base.show(stream::IO, result::MFDFAResult)
    print(
        stream,
        "MFDFAResult(q=",
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
    scaling_exponent(result::MFDFAResult) -> AbstractFloat

Generalized Hurst exponent at ``q = 2``, the standard DFA scaling exponent.

# Arguments

- `result::MFDFAResult`: a multifractal result whose moment orders include
  ``q = 2``.

# Returns

- the exponent ``h(2)``, in the result's value type.

# Throws

- `ArgumentError`: if ``q = 2`` is not among the moment orders.
"""
function scaling_exponent(result::MFDFAResult)
    index = findfirst(order_q -> isapprox(order_q, 2), result.q_values)
    index === nothing &&
        throw(ArgumentError("scaling_exponent requires q = 2 among the q values"))
    return result.generalized_hurst[index]
end

@doc doc"""
    mfdfa(series; kwargs...) -> MFDFAResult

Multifractal detrended fluctuation analysis of a one-dimensional time `series`.

The detrended segment variances are combined into q-order fluctuations for every
order in `q_values`; the log-log slope of each gives a generalized Hurst exponent
``h(q)``, from which the mass exponents ``\tau(q) = q\,h(q) - 1`` and, by a
Legendre transform, the singularity strengths ``\alpha`` and spectrum
``f(\alpha)`` are derived. At ``q = 2`` the analysis reduces to [`dfa`](@ref).

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

- `ArgumentError`: if `series` has fewer than 8 points, if fewer than two distinct
  `q` values are given, or if a scale yields a zero-variance segment (for which
  the multifractal moments are undefined).
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
    sorted_q = sort(unique(float.(q_values)))
    length(sorted_q) >= 2 ||
        throw(ArgumentError("need at least two distinct q values for MFDFA"))
    scales = Int.(collect(scales))
    profile = integrated_profile(series; demean = demean)
    fluctuations = mfdfa_fluctuations(
        profile, scales, sorted_q, detrender; overlap = overlap, bidirectional = bidirectional
    )
    # Carry the fluctuation value type through every derived quantity.
    q_values_typed = eltype(fluctuations).(sorted_q)
    fits = fit_generalized_hurst(scales, fluctuations, q_values_typed; fitrange = fitrange)
    hurst_values = [fit.exponent for fit in fits]
    mass_exponent_values = compute_mass_exponents(q_values_typed, hurst_values)
    strengths, spectrum = compute_singularity_spectrum(q_values_typed, mass_exponent_values)

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
