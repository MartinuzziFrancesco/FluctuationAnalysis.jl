module FluctuationAnalysis

using Statistics: mean, std
using Markdown: @doc_str
using ConcreteStructs: @concrete

import LinearAlgebra

include("detrending.jl")
include("profile.jl")
include("scales.jl")
include("segmentation.jl")
include("fluctuation.jl")
include("fit.jl")
include("results.jl")
include("dfa.jl")
include("mfdfa.jl")
include("dcca.jl")
include("hurst.jl")
include("moving_average.jl")
include("dma.jl")
include("mfdma.jl")

export AbstractDetrender, PolynomialDetrender, detrend
export integrated_profile, logarithmic_scales
export LogLogFit, loglog_fit
export AbstractFluctuationResult, scaling_exponent, DFAResult, dfa
export MFDFAResult, mfdfa
export DCCAResult, dcca
export AbstractHurstEstimator, DetrendedFluctuationHurst, RescaledRangeHurst
export HurstResult, hurst, hurst_exponent
export MovingAverage, DMAResult, dma
export MFDMAResult, mfdma

end # module
