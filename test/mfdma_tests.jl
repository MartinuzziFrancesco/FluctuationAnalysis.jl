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
        reference = dma(series; scales = scales)
        @test result.fluctuations[:, 1] == reference.fluctuations
        @test scaling_exponent(result) == reference.fit.exponent
    end

    @testset "theta and moving_average keywords agree" begin
        series = randn(MersenneTwister(2), 3000)
        q_values = collect(-2.0:1.0:2.0)
        @test mfdma(series; q_values = q_values, theta = 0.5).fluctuations ==
            mfdma(
            series; q_values = q_values, moving_average = MovingAverage(0.5)
        ).fluctuations
    end

    @testset "scaling_exponent requires q = 2" begin
        result = mfdma(randn(MersenneTwister(4), 2000); q_values = [-1.0, 0.0, 1.0])
        @test_throws ArgumentError scaling_exponent(result)
    end

    @testset "validation" begin
        @test_throws ArgumentError mfdma([1.0, 2.0, 3.0])
        @test_throws ArgumentError mfdma(randn(2000); q_values = [1.0])
        @test_throws ArgumentError mfdma(randn(2000); q_values = [1.0, 1.0])
        @test_throws ArgumentError mfdma(ones(100); q_values = [2.0, 4.0], scales = [8, 16])
        @test_throws ArgumentError mfdma(
            randn(100); q_values = [2.0, Inf], scales = [8, 16]
        )
        @test_throws ArgumentError mfdma(
            randn(100); q_values = [2.0, NaN], scales = [8, 16]
        )
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
        scales = logarithmic_scales(
            length(cascade); minimum_scale = 16, maximum_scale = 1024
        )
        result = mfdma(cascade; q_values = collect(-4.0:0.5:4.0), scales = scales)
        @test maximum(result.generalized_hurst) - minimum(result.generalized_hurst) > 0.3
    end
end

@testitem "mfdma: backward reproduces the binomial spectrum (Gu & Zhou 2010)" begin
    using FluctuationAnalysis

    # Deterministic two-scale binomial measure with weight p; its generalized
    # Hurst exponent is known in closed form, h(q) = (1 - log2(p^q + (1-p)^q)) / q.
    # Backward MFDMA on the raw cumulative sum (the default, eq. (1)) recovers it;
    # demeaning the profile first would bias the one-sided exponents badly, so this
    # guards that the default stays demean = false.
    p = 0.3
    measure = [1.0]
    for _ in 1:13
        refined = Vector{Float64}(undef, 2 * length(measure))
        for (index, value) in pairs(measure)
            refined[2index - 1] = value * p
            refined[2index] = value * (1 - p)
        end
        measure = refined
    end

    analytic(q) = q == 0 ? -0.5 * (log2(p) + log2(1 - p)) : (1 - log2(p^q + (1 - p)^q)) / q
    q_values = [-3.0, -1.0, 1.0, 2.0, 3.0]
    scales = [16, 32, 64, 128, 256, 512, 1024]
    result = mfdma(measure; q_values = q_values, theta = 0.0, scales = scales)
    for (index, q) in enumerate(result.q_values)
        @test isapprox(result.generalized_hurst[index], analytic(q); atol = 0.04)
    end
end
