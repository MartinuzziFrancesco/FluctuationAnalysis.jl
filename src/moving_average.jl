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

- `theta::Float64`: the window position parameter.

# Throws

- `ArgumentError`: if `theta` is outside `[0, 1]`.
"""
struct MovingAverage
    theta::Float64
    function MovingAverage(theta::Real)
        0 <= theta <= 1 ||
            throw(ArgumentError("moving-average position theta must lie in [0, 1]"))
        return new(Float64(theta))
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

"""
    window_offsets(window, theta)

Number of past and future points spanned by a moving-average window of size
`window` at position `theta`, returned as `(past, future)`.

Following Gu & Zhou (2010), the window holds `future = floor((window - 1) * theta)`
future points and `past = (window - 1) - future` past points, so that
`past + future + 1 == window` exactly. Computing `past` from `future` avoids any
floating-point disagreement between the floor and ceiling forms.
"""
function window_offsets(window::Integer, theta::Real)
    window >= 2 || throw(ArgumentError("moving-average window must be at least 2"))
    0 <= theta <= 1 || throw(ArgumentError("moving-average position theta must lie in [0, 1]"))
    future = floor(Int, (window - 1) * theta)
    past = (window - 1) - future
    return past, future
end

"""
    moving_average(profile, window, moving_average_spec)

Moving-average trend of `profile` for a window of size `window`, evaluated at
every index where the full window fits.

Returns a vector of length `length(profile) - window + 1`; entry `j` is the mean
of `profile` over the window centred (per `theta`) on index `past + j`.
"""
function moving_average(
        profile::AbstractVector{<:Real}, window::Integer, moving_average_spec::MovingAverage
    )
    profile_length = length(profile)
    window <= profile_length ||
        throw(ArgumentError("window must not exceed the profile length"))
    past, future = window_offsets(window, moving_average_spec.theta)
    valid_length = profile_length - window + 1
    trend = similar(profile, float(eltype(profile)), valid_length)
    for (output_index, center) in enumerate((past + 1):(profile_length - future))
        trend[output_index] = mean(@view profile[(center - past):(center + future)])
    end
    return trend
end

"""
    moving_average_residual(profile, window, moving_average_spec)

Residual of `profile` after subtracting its [`moving_average`](@ref) trend,
restricted to the indices where the full window fits.

Returns a vector of length `length(profile) - window + 1`, the sequence
``ε(i) = y(i) - ỹ(i)`` of Gu & Zhou (2010).
"""
function moving_average_residual(
        profile::AbstractVector{<:Real}, window::Integer, moving_average_spec::MovingAverage
    )
    past, future = window_offsets(window, moving_average_spec.theta)
    trend = moving_average(profile, window, moving_average_spec)
    valid = @view profile[(past + 1):(length(profile) - future)]
    return valid .- trend
end

"""
    moving_average_variances(profile, window, moving_average_spec)

Mean squared residual of every disjoint segment of size `window` taken from the
moving-average residual of `profile`.

The residual is partitioned into the non-overlapping segments that fit within it,
in order, and each segment's mean square ``F_v^2(n) = (1/n) Σ ε_v(i)^2`` is
returned. This mirrors steps 3–4 of the Gu & Zhou (2010) MFDMA algorithm.
"""
function moving_average_variances(
        profile::AbstractVector{<:Real}, window::Integer, moving_average_spec::MovingAverage
    )
    residual = moving_average_residual(profile, window, moving_average_spec)
    length(residual) >= window ||
        throw(ArgumentError("window $window is too large to form a residual segment"))
    segments = segment_views(residual, window; overlap = false, bidirectional = false)
    return map(segment -> mean(abs2, segment), segments)
end

"""
    dma_fluctuation_curve(profile, scales, moving_average_spec)

Root-mean-square detrending-moving-average fluctuation of `profile` at each scale
in `scales`, computed as the square root of the mean segment variance.
"""
function dma_fluctuation_curve(
        profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        moving_average_spec::MovingAverage,
    )
    fluctuations = similar(profile, float(eltype(profile)), length(scales))
    for (index, scale) in pairs(scales)
        variances = moving_average_variances(profile, scale, moving_average_spec)
        fluctuations[index] = sqrt(mean(variances))
    end
    return fluctuations
end

"""
    mfdma_fluctuations(profile, scales, q_values, moving_average_spec)

Matrix of q-order moving-average fluctuations whose entry `[scale_index, q_index]`
is the fluctuation at `scales[scale_index]` for `q_values[q_index]`, following
step 5 of Gu & Zhou (2010). Reuses [`q_order_fluctuation`](@ref).
"""
function mfdma_fluctuations(
        profile::AbstractVector{<:Real},
        scales::AbstractVector{<:Integer},
        q_values::AbstractVector{<:Real},
        moving_average_spec::MovingAverage,
    )
    value_type = float(eltype(profile))
    fluctuations = Matrix{value_type}(undef, length(scales), length(q_values))
    for (scale_index, scale) in pairs(scales)
        variances = moving_average_variances(profile, scale, moving_average_spec)
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
