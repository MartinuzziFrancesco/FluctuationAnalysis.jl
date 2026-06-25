@testitem "input types: array containers" begin
    using FluctuationAnalysis
    using Random
    using OffsetArrays

    # A correlated pair so that every method (including dcca) is well defined.
    rng = MersenneTwister(1)
    common = cumsum(randn(rng, 5000))
    a = common .+ randn(rng, 5000)
    b = common .+ randn(rng, 5000)
    scales = [8, 16, 32, 64, 128, 256]
    q = [2.0, 3.0]

    # Scaling exponents from every single-series method, as one comparable bundle.
    function exponents(series)
        return (
            dfa = scaling_exponent(dfa(series; scales = scales)),
            dma = scaling_exponent(dma(series; scales = scales)),
            mfdfa = scaling_exponent(mfdfa(series; q_values = q, scales = scales)),
            mfdma = scaling_exponent(mfdma(series; q_values = q, scales = scales)),
            hurst_dfa = hurst_exponent(hurst(series; scales = scales)),
            hurst_rs = hurst_exponent(hurst(series, RescaledRangeHurst(); scales = scales)),
        )
    end
    cross(x, y) = scaling_exponent(dcca(x, y; scales = scales))

    baseline = exponents(a)
    baseline_cross = cross(a, b)

    # Containers that hold exactly the same Float64 values as `a`/`b`.
    wrappers = (
        "view" => v -> view(v, eachindex(v)),
        "OffsetArray (0-based)" => v -> OffsetArray(v, 0:(length(v) - 1)),
        "OffsetArray (negative)" => v -> OffsetArray(v, -3:(length(v) - 4)),
    )

    @testset "$name preserves every exponent" for (name, wrap) in wrappers
        wrapped = exponents(wrap(a))
        for key in keys(baseline)
            @test isapprox(wrapped[key], baseline[key]; atol = 1.0e-10)
        end
        @test isapprox(cross(wrap(a), wrap(b)), baseline_cross; atol = 1.0e-10)
    end
end

@testitem "input types: element types" begin
    using FluctuationAnalysis
    using Random

    rng = MersenneTwister(1)
    a = cumsum(randn(rng, 5000))
    scales = [8, 16, 32, 64, 128, 256]
    baseline = scaling_exponent(dfa(a; scales = scales))

    @testset "Float32 input is accepted and close" begin
        @test isapprox(scaling_exponent(dfa(Float32.(a); scales = scales)), baseline; atol = 1.0e-2)
    end

    @testset "BigFloat input is accepted and close" begin
        @test isapprox(scaling_exponent(dfa(big.(a); scales = scales)), baseline; atol = 1.0e-6)
    end

    @testset "integer input runs and gives a finite exponent" begin
        integers = rand(MersenneTwister(2), -10:10, 5000)
        for method in (
                s -> scaling_exponent(dfa(s; scales = scales)),
                s -> scaling_exponent(dma(s; scales = scales)),
                s -> hurst_exponent(hurst(s; scales = scales)),
                s -> hurst_exponent(hurst(s, RescaledRangeHurst(); scales = scales)),
            )
            @test isfinite(method(integers))
        end
    end

    @testset "rational input runs and gives a finite exponent" begin
        rationals = (1 // 10) .* rand(MersenneTwister(3), -10:10, 5000)
        @test isfinite(scaling_exponent(dfa(rationals; scales = scales)))
    end

    @testset "element type is preserved without conversion to Float64" for T in
        (Float32, BigFloat)
        series = T.(a)
        @test eltype(dfa(series; scales = scales).fluctuations) == T
        @test scaling_exponent(dfa(series; scales = scales)) isa T
        @test eltype(dma(series; scales = scales).fluctuations) == T
        @test eltype(hurst(series; scales = scales).statistic) == T
        @test eltype(hurst(series, RescaledRangeHurst(); scales = scales).statistic) == T
        @test hurst_exponent(hurst(series; scales = scales)) isa T
    end

    @testset "BigFloat flows through the multifractal and cross methods" begin
        series = big.(a)
        other = big.(cumsum(randn(MersenneTwister(5), length(a))) .+ a)
        @test eltype(mfdfa(series; q_values = [2.0, 3.0], scales = scales).generalized_hurst) ==
            BigFloat
        @test eltype(mfdma(series; q_values = [2.0, 3.0], scales = scales).generalized_hurst) ==
            BigFloat
        @test eltype(dcca(series, other; scales = scales).correlation) == BigFloat
    end

    @testset "multifractal preserves the series type regardless of q" begin
        # Float32 data keeps Float32 even though q_values are Float64 (cast to data type).
        @test eltype(mfdfa(Float32.(a); q_values = [2.0, 3.0], scales = scales).generalized_hurst) ==
            Float32
        @test eltype(mfdma(Float32.(a); q_values = [2.0, 3.0], scales = scales).generalized_hurst) ==
            Float32
        @test eltype(mfdfa(Float32.(a); q_values = [2.0, 3.0], scales = scales).q_values) == Float32
    end
end

@testitem "input types: scale containers" begin
    using FluctuationAnalysis
    using Random
    using OffsetArrays

    a = cumsum(randn(MersenneTwister(1), 5000))
    scale_values = [8, 16, 32, 64, 128, 256]
    baseline = scaling_exponent(dfa(a; scales = scale_values))

    @testset "offset scales give the same exponent" begin
        offset_scales = OffsetArray(scale_values, 0:(length(scale_values) - 1))
        @test isapprox(scaling_exponent(dfa(a; scales = offset_scales)), baseline; atol = 1.0e-12)
    end

    @testset "stepped integer input is accepted as scales" begin
        @test isfinite(scaling_exponent(dfa(a; scales = 16:16:512)))
    end

    @testset "non-Int integer scales are accepted" begin
        @test isapprox(
            scaling_exponent(dfa(a; scales = Int32.(scale_values))), baseline; atol = 1.0e-12
        )
    end
end
