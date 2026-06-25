@testitem "mfdfa: internals" begin
    using FluctuationAnalysis
    q_order_fluctuation = FluctuationAnalysis.q_order_fluctuation

    @testset "q-order fluctuation of equal variances" begin
        variances = [4.0, 4.0, 4.0]
        for order_q in (-2.0, 0.0, 1.0, 3.0)
            @test isapprox(q_order_fluctuation(variances, order_q), 2.0; atol = 1.0e-12)
        end
    end

    @testset "near-zero q uses the logarithmic limit" begin
        variances = [1.0, 4.0, 9.0]
        reference = q_order_fluctuation(variances, 0.0)
        @test isapprox(q_order_fluctuation(variances, 1.0e-16), reference; atol = 1.0e-12)
        @test isapprox(q_order_fluctuation(variances, -1.0e-16), reference; atol = 1.0e-12)
    end

    @testset "mass exponents follow q h(q) - 1" begin
        q_values = [-1.0, 0.0, 2.0]
        hurst_values = [0.6, 0.5, 0.4]
        @test FluctuationAnalysis.compute_mass_exponents(q_values, hurst_values) ==
            q_values .* hurst_values .- 1
    end

    @testset "central difference of a line is its slope" begin
        nodes = [1.0, 2.0, 3.0, 4.0]
        @test FluctuationAnalysis.central_difference(nodes, 3.0 .* nodes .+ 1.0) == fill(3.0, 4)
    end
end

@testitem "mfdfa: result and reductions" begin
    using FluctuationAnalysis
    using Random

    @testset "result structure" begin
        q_values = collect(-3.0:1.0:3.0)
        result = mfdfa(randn(MersenneTwister(1), 4000); q_values = q_values)
        @test result isa MFDFAResult
        @test result.q_values == sort(q_values)
        @test size(result.fluctuations) == (length(result.scales), length(result.q_values))
        @test length(result.generalized_hurst) == length(result.q_values)
        @test length(result.mass_exponents) == length(result.q_values)
        @test length(result.singularity_strengths) == length(result.q_values)
        @test length(result.singularity_spectrum) == length(result.q_values)
    end

    @testset "q values are sorted and deduplicated" begin
        result = mfdfa(randn(MersenneTwister(2), 2000); q_values = [2.0, -1.0, 2.0, 0.0])
        @test result.q_values == [-1.0, 0.0, 2.0]
    end

    @testset "h(2) equals the DFA exponent" begin
        series = randn(MersenneTwister(99), 5000)
        scales = logarithmic_scales(5000; minimum_scale = 8, maximum_scale = 500)
        result = mfdfa(series; q_values = [2.0, 4.0], scales = scales)
        @test isapprox(scaling_exponent(result), dfa(series; scales = scales).fit.exponent; atol = 1.0e-10)
    end

    @testset "mass exponents are monotonically increasing in q" begin
        result = mfdfa(randn(MersenneTwister(3), 8000); q_values = collect(-3.0:1.0:3.0))
        @test issorted(result.mass_exponents)
    end

    @testset "scaling_exponent requires q = 2" begin
        result = mfdfa(randn(MersenneTwister(4), 2000); q_values = [-1.0, 0.0, 1.0])
        @test_throws ArgumentError scaling_exponent(result)
    end

    @testset "validation" begin
        @test_throws ArgumentError mfdfa([1.0, 2.0, 3.0])
        @test_throws ArgumentError mfdfa(randn(2000); q_values = [1.0])
        @test_throws ArgumentError mfdfa(randn(2000); q_values = [1.0, 1.0])
    end
end

@testitem "mfdfa: multifractality" begin
    using FluctuationAnalysis
    using Random

    @testset "white noise is approximately monofractal" begin
        series = randn(MersenneTwister(123), 30_000)
        scales = logarithmic_scales(30_000; minimum_scale = 16, maximum_scale = 3000)
        result = mfdfa(series; q_values = collect(-3.0:1.0:3.0), scales = scales)
        @test isapprox(scaling_exponent(result), 0.5; atol = 0.05)
        @test maximum(result.generalized_hurst) - minimum(result.generalized_hurst) < 0.1
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
        result = mfdfa(cascade; q_values = collect(-4.0:0.5:4.0), scales = scales)
        @test maximum(result.generalized_hurst) - minimum(result.generalized_hurst) > 0.3
    end
end
