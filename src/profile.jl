"""
    integrated_profile(series; demean = true) -> Vector{<:AbstractFloat}

Construct the integrated profile of a time `series`, the cumulative sum of the
mean-removed values. This is the first step shared by every method in the
package. The element type is `float(eltype(series))`, so a `Float32` series gives
a `Float32` profile and a `BigFloat` series a `BigFloat` profile.

# Arguments

- `series::AbstractVector{<:Real}`: the input time series; must be nonempty.

# Keywords

- `demean::Bool = true`: subtract the series mean before integrating. When
  `false` the cumulative sum is taken directly on the input values.

# Returns

- `Vector{<:AbstractFloat}`: the integrated profile, of the same length as
  `series` and element type `float(eltype(series))`.

# Throws

- `ArgumentError`: if `series` is empty.
"""
function integrated_profile(series::AbstractVector{<:Real}; demean::Bool = true)
    isempty(series) && throw(ArgumentError("input series must be nonempty"))
    values = float.(collect(series))
    centered = demean ? values .- mean(values) : values
    return cumsum(centered)
end
