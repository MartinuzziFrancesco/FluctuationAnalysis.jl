function __segment_variance(detrender::AbstractDetrender, segment::AbstractVector{<:Real})
    residuals = detrend(detrender, segment)
    return mean(abs2, residuals)
end

function __segment_covariance(
        detrender::AbstractDetrender,
        first_segment::AbstractVector{<:Real},
        second_segment::AbstractVector{<:Real},
    )
    first_residuals = detrend(detrender, first_segment)
    second_residuals = detrend(detrender, second_segment)
    return mean(first_residuals .* second_residuals)
end

function __segment_variances(
        profile::AbstractVector{<:Real},
        scale::Integer,
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    segments = __segment_views(profile, scale; overlap = overlap, bidirectional = bidirectional)
    return map(Base.Fix1(__segment_variance, detrender), segments)
end

function __segment_covariances(
        first_profile::AbstractVector{<:Real},
        second_profile::AbstractVector{<:Real},
        scale::Integer,
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    first_segments = __segment_views(
        first_profile, scale; overlap = overlap, bidirectional = bidirectional
    )
    second_segments = __segment_views(
        second_profile, scale; overlap = overlap, bidirectional = bidirectional
    )
    value_type = promote_type(float(eltype(first_profile)), float(eltype(second_profile)))
    covariances = zeros(value_type, length(first_segments))
    for index in eachindex(first_segments, second_segments)
        covariances[index] = __segment_covariance(
            detrender, first_segments[index], second_segments[index]
        )
    end
    return covariances
end

function __fluctuation_at_scale(
        profile::AbstractVector{<:Real},
        scale::Integer,
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    variances = __segment_variances(
        profile, scale, detrender; overlap = overlap, bidirectional = bidirectional
    )
    return sqrt(mean(variances))
end

function __fluctuation_curve(
        profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        detrender::AbstractDetrender;
        overlap::Bool = false,
        bidirectional::Bool = true,
    )
    smallest_allowed = __minimum_segment_length(detrender)
    fluctuations = similar(profile, float(eltype(profile)), length(scales))
    for (index, scale) in pairs(scales)
        scale >= smallest_allowed ||
            throw(ArgumentError("scale $scale is too small for the chosen detrender"))
        fluctuations[index] = __fluctuation_at_scale(
            profile, scale, detrender; overlap = overlap, bidirectional = bidirectional
        )
    end
    return fluctuations
end
