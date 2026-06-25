@testitem "Code quality (Aqua.jl)" begin
    using Aqua
    using FluctuationAnalysis
    Aqua.test_all(FluctuationAnalysis; ambiguities = false)
end

@testitem "Code linting (JET.jl)" begin
    using JET
    using FluctuationAnalysis
    JET.test_package(FluctuationAnalysis; target_modules = (FluctuationAnalysis,))
end
