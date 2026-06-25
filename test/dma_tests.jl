@testitem "dma: structure and options" begin
    using FluctuationAnalysis
    using Random

    @testset "result structure" begin
        result = dma(randn(MersenneTwister(1), 4000))
        @test result isa DMAResult
        @test length(result.scales) == length(result.fluctuations)
        @test result.moving_average isa MovingAverage
        @test scaling_exponent(result) == result.fit.exponent
    end

    @testset "theta and moving_average keywords agree" begin
        series = randn(MersenneTwister(2), 3000)
        @test dma(series; theta = 0.5).fluctuations ==
            dma(series; moving_average = MovingAverage(0.5)).fluctuations
    end

    @testset "fitrange restricts the fitted scales" begin
        series = randn(MersenneTwister(3), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 1000)
        result = dma(series; scales = scales, fitrange = (16, 256))
        @test result.fit.fitted_scale_range[1] >= 16
        @test result.fit.fitted_scale_range[2] <= 256
    end

    @testset "validation" begin
        @test_throws ArgumentError dma([1.0, 2.0, 3.0])
    end
end

@testitem "dma: scaling convergence" begin
    using FluctuationAnalysis
    using Random

    @testset "white noise gives H ≈ 0.5 for every window position" begin
        for theta in (0.0, 0.5, 1.0)
            series = randn(MersenneTwister(42), 30_000)
            scales = logarithmic_scales(30_000; minimum_scale = 16, maximum_scale = 3000)
            @test isapprox(scaling_exponent(dma(series; theta = theta, scales = scales)), 0.5; atol = 0.05)
        end
    end

    @testset "correlated signal raises the exponent above white noise" begin
        walk = cumsum(randn(MersenneTwister(5), 20_000))
        scales = logarithmic_scales(20_000; minimum_scale = 16, maximum_scale = 2000)
        @test scaling_exponent(dma(walk; theta = 0.5, scales = scales)) > 1.0
    end
end
