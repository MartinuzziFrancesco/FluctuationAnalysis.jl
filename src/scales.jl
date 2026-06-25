"""
    logarithmic_scales(series_length; kwargs...) -> Vector{Int}

Generate logarithmically spaced, distinct integer window sizes (scales) for a
series of length `series_length`. Duplicates introduced by rounding are removed,
so the result may contain fewer than `scale_count` entries. This is the default
choice of scales for every method in the package.

# Arguments

- `series_length::Integer`: length of the series to be analysed; must be at
  least 8.

# Keywords

- `minimum_scale::Integer = 4`: smallest scale; must be at least 2.
- `maximum_scale::Integer = max(minimum_scale + 1, series_length ÷ 4)`: largest
  scale; must satisfy `minimum_scale < maximum_scale ≤ series_length`.
- `scale_count::Integer = 20`: number of scales requested before deduplication;
  must be at least 2.
- `base::Real = 2.0`: base of the logarithmic spacing; must be greater than 1.

# Returns

- `Vector{Int}`: the distinct scales in increasing order.

# Throws

- `ArgumentError`: if any argument violates the bounds above, or if fewer than
  two distinct scales can be generated.
"""
function logarithmic_scales(
        series_length::Integer;
        minimum_scale::Integer = 4,
        maximum_scale::Integer = max(minimum_scale + 1, series_length ÷ 4),
        scale_count::Integer = 20,
        base::Real = 2.0,
    )
    series_length >= 8 || throw(ArgumentError("series is too short for scaling analysis"))
    minimum_scale >= 2 || throw(ArgumentError("minimum_scale must be at least 2"))
    maximum_scale > minimum_scale ||
        throw(ArgumentError("maximum_scale must exceed minimum_scale"))
    maximum_scale <= series_length ||
        throw(ArgumentError("maximum_scale must not exceed series length"))
    scale_count >= 2 || throw(ArgumentError("scale_count must be at least 2"))
    base > 1 || throw(ArgumentError("base must be greater than 1"))

    lower_exponent = log(minimum_scale) / log(base)
    upper_exponent = log(maximum_scale) / log(base)
    exponents = range(lower_exponent, upper_exponent; length = scale_count)
    scales = unique(round.(Int, base .^ exponents))
    scales = filter(scale -> minimum_scale <= scale <= maximum_scale, scales)
    length(scales) >= 2 || throw(ArgumentError("not enough distinct scales generated"))
    return scales
end
