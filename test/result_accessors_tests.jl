@testitem "result accessors: common interface" begin
    using FluctuationAnalysis
    using Random

    series = randn(MersenneTwister(51), 1024)
    scales_input = [8, 16, 32, 64, 128]
    q_input = [2.0, 3.0]
    results = (
        dfa(series; scales = scales_input),
        mfdfa(series; q_values = q_input, scales = scales_input),
        dcca(series, series; scales = scales_input),
        dma(series; scales = scales_input),
        mfdma(series; q_values = q_input, scales = scales_input),
        hurst(series; scales = scales_input),
        hurst(series, RescaledRangeHurst(); scales = scales_input),
    )

    for result in results
        @test analysis_scales(result) === result.scales
        @test fit_results(result) isa Tuple
        @test all(fit_result -> fit_result isa LogLogFit, fit_results(result))
    end

    for result in results[[1, 3, 4, 6, 7]]
        @test length(fit_results(result)) == 1
        @test only(fit_results(result)) === result.fit
    end
    for result in results[[2, 5]]
        @test fit_results(result) == Tuple(result.fits)
    end

    @test fluctuation_values(results[1]) === results[1].fluctuations
    @test fluctuation_values(results[2]) === results[2].fluctuations
    @test fluctuation_values(results[3]) === results[3].cross_fluctuations
    @test fluctuation_values(results[4]) === results[4].fluctuations
    @test fluctuation_values(results[5]) === results[5].fluctuations

    @test analysis_method(results[1]) === results[1].detrender
    @test analysis_method(results[2]) === results[2].detrender
    @test analysis_method(results[3]) === results[3].detrender
    @test analysis_method(results[4]) === results[4].moving_average
    @test analysis_method(results[5]) === results[5].moving_average
    @test analysis_method(results[6]) === results[6].estimator
    @test analysis_method(results[7]) === results[7].estimator
end

@testitem "result accessors: multifractal interface" begin
    using FluctuationAnalysis
    using Random

    series = randn(MersenneTwister(52), 1024)
    scales_input = [8, 16, 32, 64, 128]
    q_input = [-1.0, 0.0, 2.0, 3.0]

    for result in (
            mfdfa(series; q_values = q_input, scales = scales_input),
            mfdma(series; q_values = q_input, scales = scales_input),
        )
        @test moment_orders(result) === result.q_values
        @test generalized_hurst(result) === result.generalized_hurst
        @test mass_exponents(result) === result.mass_exponents
        @test singularity_strengths(result) === result.singularity_strengths
        @test singularity_spectrum(result) === result.singularity_spectrum
        @test length(moment_orders(result)) == length(fit_results(result))
        @test size(fluctuation_values(result)) ==
            (length(analysis_scales(result)), length(moment_orders(result)))
    end
end

@testitem "result accessors: DCCA and Hurst interfaces" begin
    using FluctuationAnalysis
    using Random

    random = MersenneTwister(53)
    series = randn(random, 1024)
    other_series = series .+ 0.2 .* randn(random, length(series))
    scales_input = [8, 16, 32, 64, 128]

    cross_result = dcca(series, other_series; scales = scales_input)
    @test dcca_covariances(cross_result) === cross_result.covariances
    @test dcca_correlation(cross_result) === cross_result.correlation
    marginal_fluctuations = dcca_marginal_fluctuations(cross_result)
    @test marginal_fluctuations.first === cross_result.first_fluctuations
    @test marginal_fluctuations.second === cross_result.second_fluctuations

    for estimator in (DetrendedFluctuationHurst(), RescaledRangeHurst())
        hurst_result = hurst(series, estimator; scales = scales_input)
        @test hurst_statistic(hurst_result) === hurst_result.statistic
    end
end

@testitem "result accessors: scalar types" begin
    using FluctuationAnalysis
    using Random

    source = randn(MersenneTwister(54), 256)
    scales_input = [8, 16, 32]

    for scalar_type in (Float32, Float64, BigFloat)
        series = scalar_type.(source)
        q_input = scalar_type[2, 3]
        dfa_result = dfa(series; scales = scales_input)
        multifractal_result = mfdfa(
            series; q_values = q_input, scales = scales_input
        )
        cross_result = dcca(series, series; scales = scales_input)
        hurst_result = hurst(series; scales = scales_input)

        @test eltype(fluctuation_values(dfa_result)) == scalar_type
        @test eltype(fluctuation_values(multifractal_result)) == scalar_type
        @test eltype(fluctuation_values(cross_result)) == scalar_type
        @test eltype(moment_orders(multifractal_result)) == scalar_type
        @test eltype(generalized_hurst(multifractal_result)) == scalar_type
        @test eltype(dcca_covariances(cross_result)) == scalar_type
        @test eltype(dcca_correlation(cross_result)) == scalar_type
        @test eltype(hurst_statistic(hurst_result)) == scalar_type
    end
end
