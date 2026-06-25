@testitem "mfdma: result and reductions" begin
    using FluctuationAnalysis
    using Random

    @testset "result structure" begin
        q_values = collect(-3.0:1.0:3.0)
        result = mfdma(randn(MersenneTwister(1), 4000); q_values = q_values)
        @test result isa MFDMAResult
        @test result.q_values == sort(q_values)
        @test size(result.fluctuations) == (length(result.scales), length(result.q_values))
        @test length(result.generalized_hurst) == length(result.q_values)
        @test length(result.mass_exponents) == length(result.q_values)
        @test length(result.singularity_strengths) == length(result.q_values)
        @test length(result.singularity_spectrum) == length(result.q_values)
        @test result.moving_average isa MovingAverage
    end

    @testset "q values are sorted and deduplicated" begin
        result = mfdma(randn(MersenneTwister(3), 2000); q_values = [2.0, -1.0, 2.0, 0.0])
        @test result.q_values == [-1.0, 0.0, 2.0]
    end

    @testset "h(2) equals the DMA exponent" begin
        series = randn(MersenneTwister(99), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 500)
        result = mfdma(series; q_values = [2.0, 4.0], scales = scales)
        @test isapprox(scaling_exponent(result), dma(series; scales = scales).fit.exponent; atol = 1.0e-10)
    end

    @testset "theta and moving_average keywords agree" begin
        series = randn(MersenneTwister(2), 3000)
        q_values = collect(-2.0:1.0:2.0)
        @test mfdma(series; q_values = q_values, theta = 0.5).fluctuations ==
            mfdma(series; q_values = q_values, moving_average = MovingAverage(0.5)).fluctuations
    end

    @testset "scaling_exponent requires q = 2" begin
        result = mfdma(randn(MersenneTwister(4), 2000); q_values = [-1.0, 0.0, 1.0])
        @test_throws ArgumentError scaling_exponent(result)
    end

    @testset "validation" begin
        @test_throws ArgumentError mfdma([1.0, 2.0, 3.0])
        @test_throws ArgumentError mfdma(randn(2000); q_values = [1.0])
        @test_throws ArgumentError mfdma(randn(2000); q_values = [1.0, 1.0])
    end
end

@testitem "mfdma: multifractality" begin
    using FluctuationAnalysis
    using Random

    @testset "white noise is approximately monofractal" begin
        series = randn(MersenneTwister(123), 30_000)
        scales = logarithmic_scales(30_000; minimum_scale = 16, maximum_scale = 3000)
        result = mfdma(series; q_values = collect(-3.0:1.0:3.0), scales = scales)
        @test isapprox(scaling_exponent(result), 0.5; atol = 0.06)
        @test maximum(result.generalized_hurst) - minimum(result.generalized_hurst) < 0.12
    end

    @testset "binomial cascade is multifractal" begin
        function binomial_cascade(levels, multiplier, rng)
            measure = [1.0]
            for _ in 1:levels
                refined = Vector{Float64}(undef, 2 * length(measure))
                for (index, value) in pairs(measure)
                    left, right = rand(rng) < 0.5 ? (multiplier, 1 - multiplier) :
                        (1 - multiplier, multiplier)
                    refined[2index - 1] = value * left
                    refined[2index] = value * right
                end
                measure = refined
            end
            return measure
        end

        cascade = binomial_cascade(14, 0.3, MersenneTwister(7))
        scales = logarithmic_scales(length(cascade); minimum_scale = 16, maximum_scale = 1024)
        result = mfdma(cascade; q_values = collect(-4.0:0.5:4.0), scales = scales)
        @test maximum(result.generalized_hurst) - minimum(result.generalized_hurst) > 0.3
    end
end
