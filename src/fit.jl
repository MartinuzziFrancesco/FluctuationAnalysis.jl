"""
    LogLogFit

Result of a straight-line fit in log-log coordinates, produced by
[`loglog_fit`](@ref) and embedded in every result object. The scaling exponent
of an analysis is its `exponent` field, conventionally read through
[`scaling_exponent`](@ref).

# Fields

- `exponent`: slope of the fit, the scaling exponent.
- `intercept`: intercept of the fit in log-log coordinates.
- `rsquared`: coefficient of determination of the fit.
- `fitted_scale_range::Tuple{Int,Int}`: smallest and largest scale included in
  the fit.

The numeric fields share the value type `float(eltype(fluctuations))`.
"""
@concrete struct LogLogFit
    exponent
    intercept
    rsquared
    fitted_scale_range::Tuple{Int, Int}
end

"""
    coefficient_of_determination(observed, predicted)

Coefficient of determination of `predicted` against `observed`.
"""
function coefficient_of_determination(observed::AbstractVector, predicted::AbstractVector)
    residual_sum = sum(abs2, observed .- predicted)
    total_sum = sum(abs2, observed .- mean(observed))
    total_sum == 0 && return one(residual_sum)
    return 1 - residual_sum / total_sum
end

"""
    scale_selection(scales, fitrange)

Boolean mask selecting the scales that lie within `fitrange`, given as a
`(lower, upper)` tuple. When `fitrange` is `nothing` every scale is selected.
"""
function scale_selection(scales::AbstractVector{<:Integer}, fitrange::Nothing)
    return trues(length(scales))
end

function scale_selection(
        scales::AbstractVector{<:Integer}, fitrange::Tuple{<:Integer, <:Integer}
    )
    lower, upper = fitrange
    return (scales .>= lower) .& (scales .<= upper)
end

"""
    loglog_fit(scales, fluctuations; fitrange = nothing) -> LogLogFit

Fit a straight line to `log.(fluctuations)` against `log.(scales)` by ordinary
least squares and return a [`LogLogFit`](@ref). The slope is the scaling
exponent of the analysis.

# Arguments

- `scales::AbstractVector{<:Integer}`: the window sizes; must be positive.
- `fluctuations::AbstractVector{<:Real}`: the fluctuation value at each scale;
  must be finite and positive over the fitted range.

# Keywords

- `fitrange::Union{Nothing,Tuple{<:Integer,<:Integer}} = nothing`: when given as
  `(lower, upper)`, only scales within those inclusive bounds are used; otherwise
  every scale is used.

# Returns

- `LogLogFit`: the slope, intercept, coefficient of determination, and
  the range of scales included in the fit.

# Throws

- `ArgumentError`: if `scales` and `fluctuations` differ in length, if fewer than
  two scales fall in the fitting range, or if any fitted scale or fluctuation is
  not positive and finite.
"""
function loglog_fit(
        scales::AbstractVector{<:Integer},
        fluctuations::AbstractVector{<:Real};
        fitrange::Union{Nothing, Tuple{<:Integer, <:Integer}} = nothing,
    )
    length(scales) == length(fluctuations) ||
        throw(ArgumentError("scales and fluctuations must have equal length"))

    scales = collect(scales)
    fluctuations = collect(fluctuations)
    selection = scale_selection(scales, fitrange)
    count(selection) >= 2 ||
        throw(ArgumentError("need at least two scales within the fitting range"))

    selected_scales = scales[selection]
    all(scale -> scale > 0, selected_scales) ||
        throw(ArgumentError("selected scales must be positive for a log-log fit"))
    selected_fluctuations = float.(fluctuations[selection])
    all(value -> isfinite(value) && value > 0, selected_fluctuations) || throw(
        ArgumentError(
            "selected fluctuations must be finite and positive for a log-log fit"
        ),
    )

    # Carry the fluctuation element type through the fit so no conversion happens.
    T = eltype(selected_fluctuations)
    log_scales = log.(T.(selected_scales))
    log_fluctuations = log.(selected_fluctuations)
    design = hcat(log_scales, ones(T, length(log_scales)))
    coefficients = design \ log_fluctuations
    rsquared = coefficient_of_determination(log_fluctuations, design * coefficients)

    return LogLogFit(
        coefficients[1],
        coefficients[2],
        rsquared,
        (minimum(selected_scales), maximum(selected_scales)),
    )
end
