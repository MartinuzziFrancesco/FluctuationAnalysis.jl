@testitem "hurst: internals and structure" begin
    using FluctuationAnalysis
    using Random

    @testset "rescaled range of a known segment" begin
        @test FluctuationAnalysis.rescaled_range([1.0, 2.0, 1.0, 2.0]) == 1.0
    end

    @testset "rescaled range of a constant segment is zero" begin
        @test FluctuationAnalysis.rescaled_range([3.0, 3.0, 3.0]) == 0.0
    end

    @testset "mean rescaled range averages over segments" begin
        series = [1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0]
        @test FluctuationAnalysis.mean_rescaled_range(series, 4) == 1.0
    end

    @testset "result structure and default estimator" begin
        result = hurst(randn(MersenneTwister(1), 2000))
        @test result isa HurstResult
        @test result.estimator isa DetrendedFluctuationHurst
        @test length(result.scales) == length(result.statistic)
        @test hurst_exponent(result) == result.fit.exponent
        @test scaling_exponent(result) == result.fit.exponent
    end

    @testset "DFA-based estimate matches dfa exponent" begin
        series = randn(MersenneTwister(11), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 500)
        from_hurst = hurst(series, DetrendedFluctuationHurst(); scales = scales)
        @test from_hurst.fit.exponent == dfa(series; scales = scales).fit.exponent
    end

    @testset "order keyword builds the detrender" begin
        @test DetrendedFluctuationHurst(; order = 2).detrender.order == 2
    end

    @testset "short series is rejected" begin
        @test_throws ArgumentError hurst([1.0, 2.0, 3.0])
    end
end

@testitem "hurst: estimates across signal types" begin
    using FluctuationAnalysis
    using Random

    @testset "white noise has Hurst near one half" begin
        series = randn(MersenneTwister(2024), 20_000)
        scales = logarithmic_scales(20_000; minimum_scale = 8, maximum_scale = 2000)
        @test isapprox(
            hurst_exponent(hurst(series, DetrendedFluctuationHurst(); scales = scales)),
            0.5;
            atol = 0.05,
        )
        @test isapprox(
            hurst_exponent(hurst(series, RescaledRangeHurst(); scales = scales)), 0.5; atol = 0.1
        )
    end

    @testset "persistent signal has Hurst above one half" begin
        series = cumsum(randn(MersenneTwister(77), 10_000))
        scales = logarithmic_scales(10_000; minimum_scale = 8, maximum_scale = 1000)
        @test hurst_exponent(hurst(series; scales = scales)) > 0.9
    end
end
