@testitem "detrending" begin
    using FluctuationAnalysis
    using Random

    @testset "constructor validation" begin
        @test PolynomialDetrender(2).order == 2
        @test PolynomialDetrender(0).order == 0
        @test_throws ArgumentError PolynomialDetrender(-1)
    end

    @testset "minimum segment length" begin
        @test FluctuationAnalysis.minimum_segment_length(PolynomialDetrender(0)) == 2
        @test FluctuationAnalysis.minimum_segment_length(PolynomialDetrender(1)) == 3
        @test FluctuationAnalysis.minimum_segment_length(PolynomialDetrender(3)) == 5
    end

    @testset "design matrix" begin
        design = FluctuationAnalysis.polynomial_design_matrix(5, 2)
        @test size(design) == (5, 3)
        @test design[:, 1] == ones(5)
        @test design[1, 2] == 0.0
        @test design[end, 2] == 1.0
    end

    @testset "order-0 detrender removes the mean" begin
        segment = [1.0, 5.0, 3.0, 7.0]   # mean 4.0
        residuals = detrend(PolynomialDetrender(0), segment)
        @test isapprox(sum(residuals), 0.0; atol = 1.0e-12)
        @test isapprox(residuals, segment .- 4.0; atol = 1.0e-12)
    end

    @testset "polynomial trends are removed exactly" begin
        positions = range(0.0, 1.0; length = 12)
        for order in 1:4
            coefficients = randn(MersenneTwister(order), order + 1)
            segment = sum(c * positions .^ (p - 1) for (p, c) in pairs(coefficients))
            residuals = detrend(PolynomialDetrender(order), segment)
            @test maximum(abs, residuals) < 1.0e-8
        end
    end

    @testset "linear detrender leaves curvature" begin
        positions = range(0.0, 1.0; length = 12)
        residuals = detrend(PolynomialDetrender(1), positions .^ 2)
        @test maximum(abs, residuals) > 1.0e-3
    end

    @testset "residuals are a one-based Vector" begin
        residuals = detrend(PolynomialDetrender(1), [3.0, 1.0, 4.0, 1.0, 5.0])
        @test residuals isa Vector{Float64}
        @test firstindex(residuals) == 1
    end

    @testset "segment shorter than order is rejected" begin
        @test_throws ArgumentError detrend(PolynomialDetrender(3), [1.0, 2.0, 3.0])
    end
end
