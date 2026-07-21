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
    scaling_exponent(result::AbstractFluctuationResult) -> Real

Scaling exponent of an analysis: the slope of the log-log fit held by `result`.

The default reads `result.fit.exponent`. Multifractal results override this to
return the generalized Hurst exponent at `q = 2`.

# Arguments

- `result::AbstractFluctuationResult`: any result returned by the package.

# Returns

- the scaling exponent, in the result's value type (the floated input type).
"""
scaling_exponent(result::AbstractFluctuationResult) = result.fit.exponent

function __show_scaling_summary(stream::IO, result::AbstractFluctuationResult)
    print(
        stream,
        nameof(typeof(result)),
        "(exponent=",
        round(result.fit.exponent; digits = 4),
        ", scales=",
        length(result.scales),
        ", rsquared=",
        round(result.fit.rsquared; digits = 4),
        ")",
    )
    return nothing
end

function __show_multifractal_summary(stream::IO, result::AbstractFluctuationResult)
    print(
        stream,
        nameof(typeof(result)),
        "(q=",
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

Base.show(stream::IO, result::DFAResult) = __show_scaling_summary(stream, result)
