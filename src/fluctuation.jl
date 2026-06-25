"""
    segment_variance(detrender, segment)

Mean squared residual of `segment` after detrending.
"""
function segment_variance(detrender::AbstractDetrender, segment::AbstractVector{<:Real})
    residuals = detrend(detrender, segment)
    return mean(abs2, residuals)
end

"""
    segment_covariance(detrender, first_segment, second_segment)

Mean product of the residuals of two segments after detrending each. With equal
segments this reduces to [`segment_variance`](@ref).
"""
function segment_covariance(
        detrender::AbstractDetrender,
        first_segment::AbstractVector{<:Real},
        second_segment::AbstractVector{<:Real},
    )
    first_residuals = detrend(detrender, first_segment)
    second_residuals = detrend(detrender, second_segment)
    return mean(first_residuals .* second_residuals)
end

"""
    segment_variances(profile, scale, detrender; overlap=false, bidirectional=true)

Detrended variance of every segment of `profile` at a single `scale`.
"""
function segment_variances(
        profile::AbstractVector{<:Real},
        scale::Integer,
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    segments = segment_views(profile, scale; overlap = overlap, bidirectional = bidirectional)
    return map(segment -> segment_variance(detrender, segment), segments)
end

"""
    segment_covariances(first_profile, second_profile, scale, detrender; overlap=false, bidirectional=true)

Detrended covariance of every aligned segment pair of two profiles at a single
`scale`.
"""
function segment_covariances(
        first_profile::AbstractVector{<:Real},
        second_profile::AbstractVector{<:Real},
        scale::Integer,
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    first_segments = segment_views(
        first_profile, scale; overlap = overlap, bidirectional = bidirectional
    )
    second_segments = segment_views(
        second_profile, scale; overlap = overlap, bidirectional = bidirectional
    )
    return map(
        (first_segment, second_segment) ->
        segment_covariance(detrender, first_segment, second_segment),
        first_segments,
        second_segments,
    )
end

"""
    fluctuation_at_scale(profile, scale, detrender; overlap=false, bidirectional=true)

Root-mean-square fluctuation of `profile` at a single `scale`, computed as the
square root of the mean detrended segment variance.
"""
function fluctuation_at_scale(
        profile::AbstractVector{<:Real},
        scale::Integer,
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    variances = segment_variances(
        profile, scale, detrender; overlap = overlap, bidirectional = bidirectional
    )
    return sqrt(mean(variances))
end

"""
    fluctuation_curve(profile, scales, detrender; overlap=false, bidirectional=true)

Fluctuation function evaluated over every entry of `scales`.
"""
function fluctuation_curve(
        profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    smallest_allowed = minimum_segment_length(detrender)
    fluctuations = similar(profile, float(eltype(profile)), length(scales))
    for (index, scale) in pairs(scales)
        scale >= smallest_allowed ||
            throw(ArgumentError("scale $scale is too small for the chosen detrender"))
        fluctuations[index] = fluctuation_at_scale(
            profile, scale, detrender; overlap = overlap, bidirectional = bidirectional
        )
    end
    return fluctuations
end
