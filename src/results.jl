"""
    AbstractFluctuationResult

Supertype for the result objects returned by the fluctuation analyses.

# Interface

A subtype that stores a [`LogLogFit`](@ref) in a field named `fit` inherits the
default [`scaling_exponent`](@ref); subtypes whose scaling exponent is derived
differently (such as [`MFDFAResult`](@ref)) override it.

Consumers should use [`analysis_scales`](@ref), [`fluctuation_values`](@ref),
[`fit_results`](@ref), and the method-specific result accessors instead of
depending on concrete field layout. Subtypes extend the accessors that apply to
their result data.

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

# Public interface

Use [`analysis_scales`](@ref), [`fluctuation_values`](@ref),
[`analysis_method`](@ref), [`fit_results`](@ref), and
[`scaling_exponent`](@ref) to inspect the result. The fluctuation values and fit
preserve `float(eltype(series))`, including `Float32` and `BigFloat`.

Concrete fields are implementation details and are not part of the stable public
interface.
"""
@concrete struct DFAResult <: AbstractFluctuationResult
    scales::Vector{Int}
    fluctuations
    detrender
    fit
end

Base.show(stream::IO, result::DFAResult) = __show_scaling_summary(stream, result)
