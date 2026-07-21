@testitem "dcca: reductions and structure" begin
    using FluctuationAnalysis
    using Random

    @testset "identical series reduces to DFA" begin
        series = randn(MersenneTwister(5), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 500)
        result = dcca(series, series; scales = scales)
        reference = dfa(series; scales = scales)
        @test result.cross_fluctuations == reference.fluctuations
        @test scaling_exponent(result) == reference.fit.exponent
    end

    @testset "identical series are perfectly correlated" begin
        result = dcca(randn(MersenneTwister(6), 4000), randn(MersenneTwister(6), 4000))
        @test all(==(1.0), result.correlation)
    end

    @testset "result structure" begin
        rng = MersenneTwister(12)
        common = cumsum(randn(rng, 3000))
        first_series = common .+ randn(rng, 3000)
        second_series = common .+ randn(rng, 3000)
        result = dcca(first_series, second_series)
        @test result isa DCCAResult
        count = length(result.scales)
        @test length(result.covariances) == count
        @test length(result.cross_fluctuations) == count
        @test length(result.first_fluctuations) == count
        @test length(result.second_fluctuations) == count
        @test length(result.correlation) == count
        @test result.detrender isa PolynomialDetrender
    end

    @testset "validation" begin
        @test_throws ArgumentError dcca(randn(2000), randn(1999))      # unequal length
        @test_throws ArgumentError dcca([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])  # too short
        @test_throws ArgumentError dcca(ones(100), collect(1.0:100.0); scales = [8, 16])
    end
end

@testitem "dcca: correlation behaviour" begin
    using FluctuationAnalysis
    using Random
    using Statistics

    @testset "correlation stays within bounds" begin
        result = dcca(randn(MersenneTwister(8), 6000), randn(MersenneTwister(81), 6000))
        @test all(value -> -1.0 <= value <= 1.0, result.correlation)
    end

    @testset "signed and zero-denominator correlation conventions" begin
        correlation = FluctuationAnalysis.__dcca_correlation(
            [-2.0, 2.0], [2.0, 2.0], [1.0, 1.0], [8, 16]
        )
        @test correlation == [-1.0, 1.0]
        @test_throws ArgumentError FluctuationAnalysis.__dcca_correlation(
            [0.0, 1.0], [0.0, 1.0], [1.0, 1.0], [8, 16]
        )
    end

    @testset "coupled series are strongly positively correlated" begin
        rng = MersenneTwister(10)
        common = cumsum(randn(rng, 8000))
        first_series = common .+ 0.5 .* randn(rng, 8000)
        second_series = common .+ 0.5 .* randn(rng, 8000)
        @test mean(dcca(first_series, second_series).correlation) > 0.5
    end

    @testset "independent series are weakly correlated" begin
        result = dcca(randn(MersenneTwister(20), 8000), randn(MersenneTwister(21), 8000))
        @test abs(mean(result.correlation)) < 0.3
    end

    @testset "fitrange restricts the cross-correlation fit" begin
        rng = MersenneTwister(30)
        common = cumsum(randn(rng, 8000))
        a = common .+ randn(rng, 8000)
        b = common .+ randn(rng, 8000)
        scales = logarithmic_scales(8000; minimum_scale = 8, maximum_scale = 1000)
        result = dcca(a, b; scales = scales, fitrange = (16, 256))
        @test result.fit.fitted_scale_range[1] >= 16
        @test result.fit.fitted_scale_range[2] <= 256
    end

    @testset "anti-correlated series have no cross-correlation exponent" begin
        # x and -x have everywhere-negative detrended covariance, so the exponent
        # is undefined and the coefficient is -1; dcca must direct users to it.
        series = cumsum(randn(MersenneTwister(40), 4000))
        @test_throws ArgumentError dcca(series, -series)
    end
end
