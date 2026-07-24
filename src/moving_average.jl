"""
    MovingAverage(theta = 0.0) -> MovingAverage

Moving-average detrending specification for the detrending moving average family
of methods ([`dma`](@ref), [`mfdma`](@ref)), following Gu & Zhou (2010).

The position parameter fixes where the moving-average window sits relative to the
point being detrended:

- `theta = 0` — *backward* moving average (window over past points only); the
  most accurate variant in Gu & Zhou's experiments and the package default.
- `theta = 0.5` — *centered* moving average (half past, half future); the only
  variant that removes a linear trend exactly.
- `theta = 1` — *forward* moving average (window over future points only).

Unlike an [`AbstractDetrender`](@ref), the moving average is computed globally
over the whole profile rather than per segment, so it is not a subtype of
`AbstractDetrender` and is consumed directly by [`dma`](@ref) and [`mfdma`](@ref).

# Arguments

- `theta::Real = 0.0`: window position; must lie in `[0, 1]`.

# Fields

- `theta::T`: the floated window position, preserving the precision of the input.

# Throws

- `ArgumentError`: if `theta` is outside `[0, 1]`.
"""
struct MovingAverage{T <: AbstractFloat}
    theta::T
    function MovingAverage(theta::Real)
        0 <= theta <= 1 ||
            throw(ArgumentError("moving-average position theta must lie in [0, 1]"))
        floated_theta = float(theta)
        return new{typeof(floated_theta)}(floated_theta)
    end
end

MovingAverage() = MovingAverage(0.0)

function Base.show(stream::IO, moving_average::MovingAverage)
    label = if moving_average.theta == 0
        "backward"
    elseif moving_average.theta == 1
        "forward"
    elseif moving_average.theta == 0.5
        "centered"
    else
        "theta=$(moving_average.theta)"
    end
    print(stream, "MovingAverage(", label, ")")
    return nothing
end

function __window_offsets(window::Integer, theta::Real)
    window >= 2 || throw(ArgumentError("moving-average window must be at least 2"))
    0 <= theta <= 1 ||
        throw(ArgumentError("moving-average position theta must lie in [0, 1]"))
    future = floor(Int, (window - 1) * theta)
    past = (window - 1) - future
    return past, future
end

function __moving_average_trend(
        profile::AbstractVector{<:Real}, window::Integer, moving_average_spec::MovingAverage
    )
    profile_length = length(profile)
    window <= profile_length ||
        throw(ArgumentError("window must not exceed the profile length"))
    past, future = __window_offsets(window, moving_average_spec.theta)
    valid_length = profile_length - window + 1
    trend = similar(profile, float(eltype(profile)), valid_length)
    for (output_index, center) in enumerate((past + 1):(profile_length - future))
        trend[output_index] = mean(@view profile[(center - past):(center + future)])
    end
    return trend
end

function __moving_average_residual(
        profile::AbstractVector{<:Real}, window::Integer, moving_average_spec::MovingAverage
    )
    past, future = __window_offsets(window, moving_average_spec.theta)
    trend = __moving_average_trend(profile, window, moving_average_spec)
    valid = @view profile[(past + 1):(length(profile) - future)]
    return valid .- trend
end

function __moving_average_variances(
        profile::AbstractVector{<:Real}, window::Integer, moving_average_spec::MovingAverage
    )
    residual = __moving_average_residual(profile, window, moving_average_spec)
    length(residual) >= window ||
        throw(ArgumentError("window $window is too large to form a residual segment"))
    segments = __segment_views(residual, window; overlap = false, bidirectional = false)
    return map(Base.Fix1(mean, abs2), segments)
end

function __dma_fluctuation_curve(
        profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        moving_average_spec::MovingAverage,
    )
    fluctuations = similar(profile, float(eltype(profile)), length(scales))
    for (index, scale) in pairs(scales)
        variances = __moving_average_variances(profile, scale, moving_average_spec)
        fluctuations[index] = sqrt(mean(variances))
    end
    return fluctuations
end

function __mfdma_fluctuations(
        profile::AbstractVector{T},
        scales::AbstractVector{<:Integer},
        q_values::AbstractVector{<:Real},
        moving_average_spec::MovingAverage,
    ) where {T <: Real}
    fluctuations = zeros(T, length(scales), length(q_values))
    for (scale_index, scale) in pairs(scales)
        variances = __moving_average_variances(profile, scale, moving_average_spec)
        __q_order_fluctuations!(
            view(fluctuations, scale_index, :), variances, q_values, scale
        )
    end
    return fluctuations
end
