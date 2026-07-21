@testitem "ForwardDiff: exported analyses differentiate series values" begin
    using FluctuationAnalysis
    using ForwardDiff
    using Random

    rng = MersenneTwister(41)
    series = randn(rng, 256)
    other_series = series .+ 0.2 .* randn(rng, 256)
    direction = randn(rng, 256)
    scales = [8, 16, 32, 64]

    analyses = (
        dfa = values -> scaling_exponent(dfa(values; scales = scales)),
        mfdfa = values -> scaling_exponent(
            mfdfa(values; q_values = [2.0, 3.0], scales = scales)
        ),
        dcca = values -> scaling_exponent(dcca(values, other_series; scales = scales)),
        dma = values -> scaling_exponent(dma(values; scales = scales)),
        mfdma = values -> scaling_exponent(
            mfdma(values; q_values = [2.0, 3.0], scales = scales)
        ),
        hurst_dfa = values -> hurst_exponent(hurst(values; scales = scales)),
        hurst_rescaled_range = values -> hurst_exponent(
            hurst(values, RescaledRangeHurst(); scales = scales)
        ),
    )

    step = 1.0e-5
    for (name, analysis) in pairs(analyses)
        forward_derivative = ForwardDiff.derivative(
            multiplier -> analysis(series .+ multiplier .* direction), 0.0
        )
        finite_difference =
            (analysis(series .+ step .* direction) - analysis(series .- step .* direction)) /
            (2step)
        @testset "$name" begin
            @test isfinite(forward_derivative)
            @test !iszero(forward_derivative)
            @test isapprox(forward_derivative, finite_difference; rtol = 2.0e-4, atol = 1.0e-7)
        end
    end
end

@testitem "ForwardDiff: shared primitives preserve derivatives" begin
    using FluctuationAnalysis
    using ForwardDiff

    series = collect(range(-1.0, 2.0; length = 32)) .^ 2
    profile_derivative = ForwardDiff.derivative(
        offset -> sum(integrated_profile(series .+ offset; demean = false)), 0.0
    )
    detrend_derivative = ForwardDiff.derivative(
        offset -> sum(abs2, detrend(PolynomialDetrender(2), series .+ offset)), 0.0
    )

    @test profile_derivative == length(series) * (length(series) + 1) / 2
    @test isfinite(detrend_derivative)
end
