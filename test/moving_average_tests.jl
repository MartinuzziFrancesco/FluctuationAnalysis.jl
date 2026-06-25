@testitem "moving_average: primitives" begin
    using FluctuationAnalysis
    using Random
    using Statistics
    FA = FluctuationAnalysis

    @testset "MovingAverage construction and validation" begin
        @test MovingAverage().theta == 0.0
        @test MovingAverage(0.5).theta == 0.5
        @test_throws ArgumentError MovingAverage(-0.1)
        @test_throws ArgumentError MovingAverage(1.1)
    end

    @testset "window_offsets split past and future" begin
        @test FA.window_offsets(5, 0.0) == (4, 0)   # backward: all past
        @test FA.window_offsets(5, 1.0) == (0, 4)   # forward: all future
        @test FA.window_offsets(5, 0.5) == (2, 2)   # centered
        @test FA.window_offsets(4, 0.5) == (2, 1)
        for n in 2:40, theta in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
            past, future = FA.window_offsets(n, theta)
            @test past + future == n - 1
            @test future == floor(Int, (n - 1) * theta)
        end
        @test_throws ArgumentError FA.window_offsets(1, 0.0)
    end

    @testset "moving average of a constant is the constant" begin
        trend = FA.moving_average(fill(3.0, 12), 4, MovingAverage(0.0))
        @test length(trend) == 12 - 4 + 1
        @test all(value -> value ≈ 3.0, trend)
    end

    @testset "moving average tracks a linear ramp with the expected lag" begin
        ramp = collect(1.0:10.0)
        @test FA.moving_average(ramp, 4, MovingAverage(0.0)) ≈ collect(4.0:10.0) .- 1.5
        @test FA.moving_average(ramp, 4, MovingAverage(1.0)) ≈ collect(1.0:7.0) .+ 1.5
    end

    @testset "centered MA removes a line; one-sided MA leaves a constant offset" begin
        ramp = collect(1.0:10.0)
        @test maximum(abs, FA.moving_average_residual(ramp, 5, MovingAverage(0.5))) < 1.0e-10
        @test all(value -> value ≈ 1.5, FA.moving_average_residual(ramp, 4, MovingAverage(0.0)))
        @test all(value -> value ≈ -1.5, FA.moving_average_residual(ramp, 4, MovingAverage(1.0)))
    end

    @testset "segment variances vanish for a perfectly detrended signal" begin
        ramp = collect(1.0:100.0)
        @test maximum(FA.moving_average_variances(ramp, 5, MovingAverage(0.5))) < 1.0e-9
        @test maximum(FA.moving_average_variances(fill(2.0, 100), 5, MovingAverage(0.0))) < 1.0e-12
        @test length(FA.moving_average_variances(ramp, 5, MovingAverage(0.0))) ==
            fld(100 - 5 + 1, 5)
    end

    @testset "q = 2 fluctuation equals the RMS fluctuation curve" begin
        profile = cumsum(randn(MersenneTwister(7), 2000))
        scales = [8, 16, 32, 64]
        curve = FA.dma_fluctuation_curve(profile, scales, MovingAverage(0.0))
        @test curve ≈ FA.mfdma_fluctuations(profile, scales, [2.0], MovingAverage(0.0))[:, 1]
    end
end

@testitem "moving_average: faithful to Gu & Zhou (2010)" begin
    using FluctuationAnalysis
    using Random
    using Statistics
    FA = FluctuationAnalysis

    # Independent reference implementation straight from eqs (1)-(6).
    function reference_fluctuation(series, window, order_q, theta)
        cumulative = cumsum(series)
        future = floor(Int, (window - 1) * theta)
        past = (window - 1) - future
        residual = Float64[]
        for index in (past + 1):(length(cumulative) - future)
            trend = mean(@view cumulative[(index - past):(index + future)])
            push!(residual, cumulative[index] - trend)
        end
        segment_count = fld(length(residual), window)
        rms = Float64[]
        for segment in 1:segment_count
            piece = @view residual[((segment - 1) * window + 1):(segment * window)]
            push!(rms, sqrt(mean(abs2, piece)))
        end
        return order_q == 0 ? exp(mean(log.(rms))) : (mean(rms .^ order_q))^(1 / order_q)
    end

    series = randn(MersenneTwister(123), 1500)
    centered = series .- mean(series)        # integrated_profile(series; demean=true) input
    profile = cumsum(centered)
    for theta in (0.0, 0.5, 1.0), window in (10, 17, 33, 64), order_q in (-3.0, 0.0, 2.0, 4.0)
        mine = FA.mfdma_fluctuations(profile, [window], [order_q], MovingAverage(theta))[1, 1]
        @test isapprox(mine, reference_fluctuation(centered, window, order_q, theta); rtol = 1.0e-12)
    end
end
