@testitem "moving_average: primitives" begin
    using FluctuationAnalysis
    using Random
    using Statistics
    const FA = FluctuationAnalysis

    @testset "MovingAverage construction and validation" begin
        @test MovingAverage().theta == 0.0
        @test MovingAverage(0.5).theta == 0.5
        @test MovingAverage(Float32(0.5)).theta isa Float32
        @test MovingAverage(big"0.5").theta isa BigFloat
        @test_throws ArgumentError MovingAverage(-0.1)
        @test_throws ArgumentError MovingAverage(1.1)
    end

    @testset "__window_offsets split past and future" begin
        @test FA.__window_offsets(5, 0.0) == (4, 0)   # backward: all past
        @test FA.__window_offsets(5, 1.0) == (0, 4)   # forward: all future
        @test FA.__window_offsets(5, 0.5) == (2, 2)   # centered
        @test FA.__window_offsets(4, 0.5) == (2, 1)
        for n in 2:40, theta in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
            past, future = FA.__window_offsets(n, theta)
            @test past + future == n - 1
            @test future == floor(Int, (n - 1) * theta)
        end
        @test_throws ArgumentError FA.__window_offsets(1, 0.0)
    end

    @testset "moving average of a constant is the constant" begin
        trend = FA.__moving_average_trend(fill(3.0, 12), 4, MovingAverage(0.0))
        @test length(trend) == 12 - 4 + 1
        @test all(value -> value ≈ 3.0, trend)
    end

    @testset "moving average tracks a linear ramp with the expected lag" begin
        ramp = collect(1.0:10.0)
        @test FA.__moving_average_trend(ramp, 4, MovingAverage(0.0)) ≈
            collect(4.0:10.0) .- 1.5
        @test FA.__moving_average_trend(ramp, 4, MovingAverage(1.0)) ≈
            collect(1.0:7.0) .+ 1.5
    end

    @testset "centered MA removes a line; one-sided MA leaves a constant offset" begin
        ramp = collect(1.0:10.0)
        centered = FA.__moving_average_residual(ramp, 5, MovingAverage(0.5))
        backward = FA.__moving_average_residual(ramp, 4, MovingAverage(0.0))
        forward = FA.__moving_average_residual(ramp, 4, MovingAverage(1.0))
        @test maximum(abs, centered) < 1.0e-10
        @test all(value -> value ≈ 1.5, backward)
        @test all(value -> value ≈ -1.5, forward)
    end

    @testset "segment variances vanish for a perfectly detrended signal" begin
        ramp = collect(1.0:100.0)
        @test maximum(FA.__moving_average_variances(ramp, 5, MovingAverage(0.5))) < 1.0e-9
        variances = FA.__moving_average_variances(
            fill(2.0, 100), 5, MovingAverage(0.0)
        )
        @test maximum(variances) < 1.0e-12
        @test length(FA.__moving_average_variances(ramp, 5, MovingAverage(0.0))) ==
            fld(100 - 5 + 1, 5)
    end

    @testset "q = 2 fluctuation equals the RMS fluctuation curve" begin
        profile = cumsum(randn(MersenneTwister(7), 2000))
        scales = [8, 16, 32, 64]
        curve = FA.__dma_fluctuation_curve(profile, scales, MovingAverage(0.0))
        multifractal = FA.__mfdma_fluctuations(
            profile, scales, [2.0], MovingAverage(0.0)
        )
        @test curve == multifractal[:, 1]
    end
end

@testitem "moving_average: faithful to Gu & Zhou (2010)" begin
    using FluctuationAnalysis
    using Random
    using Statistics
    const FA = FluctuationAnalysis

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
    theta_values = (0.0, 0.5, 1.0)
    windows = (10, 17, 33, 64)
    q_values = (-3.0, 0.0, 2.0, 4.0)
    for theta in theta_values, window in windows, order_q in q_values
        mine = FA.__mfdma_fluctuations(
            profile, [window], [order_q], MovingAverage(theta)
        )[1, 1]
        reference = reference_fluctuation(centered, window, order_q, theta)
        @test isapprox(mine, reference; rtol = 1.0e-12)
    end
end
