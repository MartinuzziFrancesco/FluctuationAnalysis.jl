"""
    AbstractFluctuationResult

Supertype for the result objects returned by the fluctuation analyses.

# Interface

A subtype that stores a [`LogLogFit`](@ref) in a field named `fit` inherits the
default [`scaling_exponent`](@ref); subtypes whose scaling exponent is derived
differently (such as [`MFDFAResult`](@ref)) override it.

See also [`DFAResult`](@ref), [`MFDFAResult`](@ref), [`DCCAResult`](@ref),
[`DMAResult`](@ref), [`MFDMAResult`](@ref), and [`HurstResult`](@ref).
"""
abstract type AbstractFluctuationResult end

"""
    scaling_exponent(result::AbstractFluctuationResult) -> AbstractFloat

Scaling exponent of an analysis: the slope of the log-log fit held by `result`.

The default reads `result.fit.exponent`. Multifractal results override this to
return the generalized Hurst exponent at `q = 2`.

# Arguments

- `result::AbstractFluctuationResult`: any result returned by the package.

# Returns

- the scaling exponent, in the result's value type (the floated input type).
"""
scaling_exponent(result::AbstractFluctuationResult) = result.fit.exponent

"""
    DFAResult

Result of a detrended fluctuation analysis, returned by [`dfa`](@ref).

# Fields

- `scales::Vector{Int}`: window sizes at which the fluctuation was evaluated.
- `fluctuations::Vector{T}`: fluctuation function values, one per scale.
- `detrender::AbstractDetrender`: detrender used to remove local trends.
- `fit::LogLogFit{T}`: the fit holding the scaling exponent and metadata.

The fluctuations and fit share the value type `float(eltype(series))`, so the
result preserves the precision of the input (`Float32`, `Float64`, `BigFloat`, ...).
"""
@concrete struct DFAResult <: AbstractFluctuationResult
    scales::Vector{Int}
    fluctuations
    detrender
    fit
end

function Base.show(stream::IO, result::DFAResult)
    print(
        stream,
        "DFAResult(exponent=",
        round(result.fit.exponent; digits = 4),
        ", scales=",
        length(result.scales),
        ", rsquared=",
        round(result.fit.rsquared; digits = 4),
        ")",
    )
    return nothing
end
