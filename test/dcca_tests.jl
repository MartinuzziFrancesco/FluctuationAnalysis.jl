@testitem "dcca: reductions and structure" begin
    using FluctuationAnalysis
    using Random

    @testset "identical series reduces to DFA" begin
        series = randn(MersenneTwister(5), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 500)
        result = dcca(series, series; scales = scales)
        reference = dfa(series; scales = scales)
        @test isapprox(result.cross_fluctuations, reference.fluctuations; atol = 1.0e-10)
        @test isapprox(scaling_exponent(result), reference.fit.exponent; atol = 1.0e-10)
    end

    @testset "identical series are perfectly correlated" begin
        result = dcca(randn(MersenneTwister(6), 4000), randn(MersenneTwister(6), 4000))
        @test all(value -> isapprox(value, 1.0; atol = 1.0e-10), result.correlation)
    end

    @testset "result structure" begin
        result = dcca(randn(MersenneTwister(12), 3000), randn(MersenneTwister(13), 3000))
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
    end
end

@testitem "dcca: correlation behaviour" begin
    using FluctuationAnalysis
    using Random
    using Statistics

    @testset "correlation stays within bounds" begin
        result = dcca(randn(MersenneTwister(8), 6000), randn(MersenneTwister(81), 6000))
        @test all(value -> -1.0 - 1.0e-8 <= value <= 1.0 + 1.0e-8, result.correlation)
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
