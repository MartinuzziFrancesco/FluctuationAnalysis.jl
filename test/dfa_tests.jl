@testitem "dfa: structure and options" begin
    using FluctuationAnalysis
    using Random

    @testset "result structure" begin
        result = dfa(randn(MersenneTwister(1), 2000))
        @test result isa DFAResult
        @test length(result.scales) == length(result.fluctuations)
        @test result.detrender isa PolynomialDetrender
        @test result.fit isa LogLogFit
        @test result.scales isa Vector{Int}
    end

    @testset "order keyword builds the detrender" begin
        @test dfa(randn(MersenneTwister(2), 2000); order = 2).detrender.order == 2
    end

    @testset "explicit detrender overrides order" begin
        result = dfa(
            randn(MersenneTwister(3), 2000);
            order = 1,
            detrender = PolynomialDetrender(3),
            scales = [8, 16, 32, 64],
        )
        @test result.detrender.order == 3
    end

    @testset "fitrange restricts the fitted scales" begin
        series = randn(MersenneTwister(4), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 1000)
        result = dfa(series; scales = scales, fitrange = (16, 256))
        @test result.fit.fitted_scale_range[1] >= 16
        @test result.fit.fitted_scale_range[2] <= 256
    end

    @testset "overlap and bidirectional are accepted" begin
        series = randn(MersenneTwister(5), 3000)
        scales = [8, 16, 32, 64, 128]
        @test dfa(series; scales = scales, overlap = true) isa DFAResult
        @test dfa(series; scales = scales, bidirectional = false) isa DFAResult
    end

    @testset "validation" begin
        @test_throws ArgumentError dfa([1.0, 2.0, 3.0])               # too short
        @test_throws ArgumentError dfa(fill(1.0, 100); scales = [8, 16, 32])  # zero fluctuation
    end
end

@testitem "dfa: scaling convergence" begin
    using FluctuationAnalysis
    using Random

    @testset "white noise scales near one half" begin
        series = randn(MersenneTwister(42), 20_000)
        scales = logarithmic_scales(20_000; minimum_scale = 8, maximum_scale = 2000)
        @test isapprox(scaling_exponent(dfa(series; scales = scales)), 0.5; atol = 0.05)
    end

    @testset "integrated white noise scales near three halves" begin
        series = cumsum(randn(MersenneTwister(7), 20_000))
        scales = logarithmic_scales(20_000; minimum_scale = 8, maximum_scale = 2000)
        @test isapprox(scaling_exponent(dfa(series; scales = scales)), 1.5; atol = 0.05)
    end

    @testset "anti-correlated signal scales below one half" begin
        # Differencing white noise produces an anti-persistent series (H < 0.5).
        noise = randn(MersenneTwister(9), 20_001)
        series = diff(noise)
        scales = logarithmic_scales(20_000; minimum_scale = 8, maximum_scale = 2000)
        @test scaling_exponent(dfa(series; scales = scales)) < 0.5
    end
end
