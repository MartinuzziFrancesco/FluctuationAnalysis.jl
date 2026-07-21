"""
    AbstractDetrender

Supertype for local trend estimators used by the fluctuation analyses.

A concrete detrender removes an estimated trend from a segment of a profile and
returns the residuals. Sharing this interface lets every per-segment method (DFA,
MFDFA, DCCA, and the Hurst estimators) reuse the same segmentation, fluctuation,
and fitting machinery.

# Interface

A concrete subtype `D <: AbstractDetrender` must implement:

- `detrend(detrender::D, segment)`: return the residuals of `segment` after
  removing the estimated trend.

See also [`PolynomialDetrender`](@ref) and [`detrend`](@ref).
"""
abstract type AbstractDetrender end

"""
    PolynomialDetrender(order) -> PolynomialDetrender

Detrender that removes a least-squares polynomial trend of degree `order` from
each segment. `order = 0` removes the mean, `order = 1` a straight line, and so
on. This is the detrender used by [`dfa`](@ref) and the other per-segment
methods.

# Arguments

- `order::Integer`: degree of the polynomial trend; must be nonnegative.

# Throws

- `ArgumentError`: if `order` is negative.
"""
struct PolynomialDetrender <: AbstractDetrender
    order::Int
    function PolynomialDetrender(order::Integer)
        order >= 0 ||
            throw(ArgumentError("polynomial detrending order must be nonnegative"))
        return new(Int(order))
    end
end

__minimum_segment_length(detrender::PolynomialDetrender) = detrender.order + 2

function __polynomial_design_matrix(
        segment_length::Integer, order::Integer, ::Type{T} = Float64
    ) where {T <: Real}
    positions = range(zero(T), one(T); length = segment_length)
    design = zeros(T, segment_length, order + 1)
    for power in 0:order
        design[:, power + 1] = positions .^ power
    end
    return design
end

"""
    detrend(detrender::AbstractDetrender, segment) -> Vector

Return the residuals of `segment` after removing the trend estimated by
`detrender`.

Part of the [`AbstractDetrender`](@ref) interface. For a [`PolynomialDetrender`](@ref)
the trend is the least-squares polynomial fit of the configured degree.

# Arguments

- `detrender::AbstractDetrender`: the detrender supplying the trend model.
- `segment::AbstractVector{<:Real}`: the contiguous piece of profile to detrend.

# Returns

- `Vector{<:Real}`: the residuals `segment .- trend`, of the same length
  and (floated) element type as `segment`.

# Throws

- `ArgumentError`: if `segment` is not longer than the polynomial degree.
"""
function detrend(detrender::PolynomialDetrender, segment::AbstractVector{<:Real})
    segment_length = length(segment)
    detrender.order < segment_length ||
        throw(ArgumentError("segment length must exceed the detrending order"))
    values = float.(collect(segment))
    design = __polynomial_design_matrix(segment_length, detrender.order, eltype(values))
    coefficients = design \ values
    return values .- design * coefficients
end
