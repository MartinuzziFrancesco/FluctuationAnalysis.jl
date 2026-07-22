@testitem "numerical reliability: additive offsets" begin
    using FluctuationAnalysis
    using Random

    random = MersenneTwister(41)
    series = randn(random, 2048)
    other_series = series .+ 0.2 .* randn(random, length(series))
    scales = [8, 16, 32, 64, 128, 256]
    q_values = [2.0, 3.0]

    shifted_series = series .+ 17.0
    shifted_other_series = other_series .- 9.0

    @test dfa(shifted_series; scales = scales).fluctuations ≈
        dfa(series; scales = scales).fluctuations
    @test mfdfa(shifted_series; q_values = q_values, scales = scales).fluctuations ≈
        mfdfa(series; q_values = q_values, scales = scales).fluctuations
    @test dma(shifted_series; scales = scales, demean = true).fluctuations ≈
        dma(series; scales = scales, demean = true).fluctuations
    @test mfdma(
        shifted_series; q_values = q_values, scales = scales, demean = true
    ).fluctuations ≈
        mfdma(series; q_values = q_values, scales = scales, demean = true).fluctuations

    reference_cross = dcca(series, other_series; scales = scales)
    shifted_cross = dcca(shifted_series, shifted_other_series; scales = scales)
    @test shifted_cross.covariances ≈ reference_cross.covariances
    @test shifted_cross.correlation ≈ reference_cross.correlation

    for estimator in (DetrendedFluctuationHurst(), RescaledRangeHurst())
        reference = hurst(series, estimator; scales = scales)
        shifted = hurst(shifted_series, estimator; scales = scales)
        @test shifted.statistic ≈ reference.statistic
        @test hurst_exponent(shifted) ≈ hurst_exponent(reference)
    end
end

@testitem "numerical reliability: signal scaling" begin
    using FluctuationAnalysis
    using Random

    random = MersenneTwister(42)
    series = randn(random, 2048)
    other_series = series .+ 0.2 .* randn(random, length(series))
    scales = [8, 16, 32, 64, 128, 256]
    q_values = [2.0, 3.0]
    multiplier = -7.0
    magnitude = abs(multiplier)

    reference_dfa = dfa(series; scales = scales)
    scaled_dfa = dfa(multiplier .* series; scales = scales)
    @test scaled_dfa.fluctuations ≈ magnitude .* reference_dfa.fluctuations
    @test scaling_exponent(scaled_dfa) ≈ scaling_exponent(reference_dfa)

    reference_mfdfa = mfdfa(series; q_values = q_values, scales = scales)
    scaled_mfdfa = mfdfa(multiplier .* series; q_values = q_values, scales = scales)
    @test scaled_mfdfa.fluctuations ≈ magnitude .* reference_mfdfa.fluctuations
    @test scaled_mfdfa.generalized_hurst ≈ reference_mfdfa.generalized_hurst

    reference_dma = dma(series; scales = scales)
    scaled_dma = dma(multiplier .* series; scales = scales)
    @test scaled_dma.fluctuations ≈ magnitude .* reference_dma.fluctuations
    @test scaling_exponent(scaled_dma) ≈ scaling_exponent(reference_dma)

    reference_mfdma = mfdma(series; q_values = q_values, scales = scales)
    scaled_mfdma = mfdma(multiplier .* series; q_values = q_values, scales = scales)
    @test scaled_mfdma.fluctuations ≈ magnitude .* reference_mfdma.fluctuations
    @test scaled_mfdma.generalized_hurst ≈ reference_mfdma.generalized_hurst

    reference_cross = dcca(series, other_series; scales = scales)
    scaled_cross = dcca(
        multiplier .* series, multiplier .* other_series; scales = scales
    )
    @test scaled_cross.covariances ≈ multiplier^2 .* reference_cross.covariances
    @test scaled_cross.cross_fluctuations ≈ magnitude .* reference_cross.cross_fluctuations
    @test scaled_cross.correlation ≈ reference_cross.correlation
    @test scaling_exponent(scaled_cross) ≈ scaling_exponent(reference_cross)

    reference_rescaled_range = hurst(series, RescaledRangeHurst(); scales = scales)
    scaled_rescaled_range = hurst(
        multiplier .* series, RescaledRangeHurst(); scales = scales
    )
    @test scaled_rescaled_range.statistic ≈ reference_rescaled_range.statistic
    @test hurst_exponent(scaled_rescaled_range) ≈ hurst_exponent(reference_rescaled_range)
end

@testitem "numerical reliability: conditioning and length" begin
    using FluctuationAnalysis
    using Random

    random = MersenneTwister(43)
    scales = [8, 16, 32, 64, 128]
    q_values = [2.0, 3.0]
    nearly_constant = 1.0 .+ 1.0e-10 .* randn(random, 1024)

    results = (
        dfa(nearly_constant; scales = scales),
        mfdfa(nearly_constant; q_values = q_values, scales = scales),
        dma(nearly_constant; scales = scales, demean = true),
        mfdma(nearly_constant; q_values = q_values, scales = scales, demean = true),
        dcca(nearly_constant, nearly_constant; scales = scales),
        hurst(nearly_constant; scales = scales),
        hurst(nearly_constant, RescaledRangeHurst(); scales = scales),
    )
    @test all(result -> isfinite(scaling_exponent(result)), results)

    boundary_series = [0.0, 1.0, -1.0, 2.0, -2.0, 1.0, 0.0, -1.0]
    boundary_scales = [3, 4]
    @test isfinite(scaling_exponent(dfa(boundary_series; scales = boundary_scales)))
    @test isfinite(scaling_exponent(dma(boundary_series; scales = boundary_scales)))
    @test isfinite(
        scaling_exponent(
            mfdfa(boundary_series; q_values = q_values, scales = boundary_scales)
        ),
    )
    @test isfinite(
        scaling_exponent(
            mfdma(boundary_series; q_values = q_values, scales = boundary_scales)
        ),
    )
    @test isfinite(
        scaling_exponent(dcca(boundary_series, boundary_series; scales = boundary_scales)),
    )
    @test isfinite(hurst_exponent(hurst(boundary_series; scales = boundary_scales)))
    @test isfinite(
        hurst_exponent(
            hurst(boundary_series, RescaledRangeHurst(); scales = boundary_scales)
        ),
    )
end

@testitem "numerical reliability: scalar types and exact reductions" begin
    using FluctuationAnalysis
    using Random

    source = randn(MersenneTwister(44), 512)
    scales = [8, 16, 32, 64]

    for scalar_type in (Float32, Float64, BigFloat)
        series = scalar_type.(source)
        q_values = scalar_type[2, 3]

        dfa_result = dfa(series; scales = scales)
        mfdfa_result = mfdfa(series; q_values = q_values, scales = scales)
        dcca_result = dcca(series, series; scales = scales)
        dma_result = dma(series; scales = scales)
        mfdma_result = mfdma(series; q_values = q_values, scales = scales)
        dfa_hurst = hurst(series; scales = scales)
        rescaled_range_hurst = hurst(series, RescaledRangeHurst(); scales = scales)

        @test eltype(dfa_result.fluctuations) == scalar_type
        @test eltype(mfdfa_result.fluctuations) == scalar_type
        @test eltype(dcca_result.covariances) == scalar_type
        @test eltype(dma_result.fluctuations) == scalar_type
        @test eltype(mfdma_result.fluctuations) == scalar_type
        @test eltype(dfa_hurst.statistic) == scalar_type
        @test eltype(rescaled_range_hurst.statistic) == scalar_type

        @test dcca_result.cross_fluctuations == dfa_result.fluctuations
        @test scaling_exponent(dcca_result) == scaling_exponent(dfa_result)
        @test mfdfa_result.fluctuations[:, 1] == dfa_result.fluctuations
        @test scaling_exponent(mfdfa_result) == scaling_exponent(dfa_result)
        @test mfdma_result.fluctuations[:, 1] == dma_result.fluctuations
        @test scaling_exponent(mfdma_result) == scaling_exponent(dma_result)
    end
end
