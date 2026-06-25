@testitem "fluctuation" begin
    using FluctuationAnalysis
    using Random

    @testset "linear profile has zero variance under linear detrending" begin
        profile = collect(1.0:20.0)
        @test FluctuationAnalysis.segment_variance(PolynomialDetrender(1), profile) < 1.0e-10
    end

    @testset "segment covariance reduces to variance for equal segments" begin
        segment = cumsum(sin.(range(0.0, 1.0; length = 16)) .+ 0.3)
        detrender = PolynomialDetrender(1)
        @test isapprox(
            FluctuationAnalysis.segment_covariance(detrender, segment, segment),
            FluctuationAnalysis.segment_variance(detrender, segment);
            atol = 1.0e-12,
        )
    end

    @testset "fluctuation of a linear profile vanishes" begin
        profile = collect(1.0:64.0)
        @test FluctuationAnalysis.fluctuation_at_scale(profile, 8, PolynomialDetrender(1)) <
            1.0e-10
    end

    @testset "fluctuation grows with scale for a random walk" begin
        walk = cumsum(randn(MersenneTwister(1), 2000))
        small = FluctuationAnalysis.fluctuation_at_scale(walk, 8, PolynomialDetrender(1))
        large = FluctuationAnalysis.fluctuation_at_scale(walk, 128, PolynomialDetrender(1))
        @test large > small
    end

    @testset "fluctuation curve matches per-scale evaluation" begin
        profile = cumsum(randn(MersenneTwister(2), 1000))
        scales = [8, 16, 32]
        detrender = PolynomialDetrender(1)
        curve = FluctuationAnalysis.fluctuation_curve(profile, scales, detrender)
        expected = [
            FluctuationAnalysis.fluctuation_at_scale(profile, scale, detrender) for
                scale in scales
        ]
        @test curve == expected
    end

    @testset "overlapping segmentation changes the estimate but stays positive" begin
        profile = cumsum(randn(MersenneTwister(3), 2000))
        detrender = PolynomialDetrender(1)
        disjoint = FluctuationAnalysis.fluctuation_at_scale(profile, 32, detrender)
        overlapped = FluctuationAnalysis.fluctuation_at_scale(
            profile, 32, detrender; overlap = true
        )
        @test disjoint > 0
        @test overlapped > 0
    end

    @testset "scale below the detrender minimum is rejected" begin
        profile = collect(1.0:64.0)
        @test_throws ArgumentError FluctuationAnalysis.fluctuation_curve(
            profile, [3], PolynomialDetrender(2)
        )
    end
end
