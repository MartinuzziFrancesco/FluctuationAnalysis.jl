@testitem "degenerate domains: exported analysis boundaries" begin
    using FluctuationAnalysis
    using Random

    series = randn(MersenneTwister(101), 100)
    other_series = series .+ 0.1 .* randn(MersenneTwister(102), 100)

    @testset "scales producing no valid segments are rejected" begin
        @test_throws ArgumentError dfa(series; scales = [8, 101])
        @test_throws ArgumentError mfdfa(series; q_values = [2.0, 3.0], scales = [8, 101])
        @test_throws ArgumentError dcca(series, other_series; scales = [8, 101])
        @test_throws ArgumentError hurst(series; scales = [8, 101])
        @test_throws ArgumentError hurst(series, RescaledRangeHurst(); scales = [8, 101])
        @test_throws ArgumentError dma(series; scales = [8, 51])
        @test_throws ArgumentError mfdma(
            series; q_values = [2.0, 3.0], scales = [8, 51]
        )
    end

    @testset "fit ranges with fewer than two scales are rejected" begin
        fitrange = (9, 15)
        scales = [8, 16]
        @test_throws ArgumentError dfa(series; scales = scales, fitrange = fitrange)
        @test_throws ArgumentError dma(series; scales = scales, fitrange = fitrange)
        @test_throws ArgumentError mfdfa(
            series; q_values = [2.0, 3.0], scales = scales, fitrange = fitrange
        )
        @test_throws ArgumentError mfdma(
            series; q_values = [2.0, 3.0], scales = scales, fitrange = fitrange
        )
        @test_throws ArgumentError hurst(series; scales = scales, fitrange = fitrange)
        @test_throws ArgumentError hurst(
            series, RescaledRangeHurst(); scales = scales, fitrange = fitrange
        )
        @test_throws ArgumentError dcca(
            series, other_series; scales = scales, fitrange = fitrange
        )
    end
end
