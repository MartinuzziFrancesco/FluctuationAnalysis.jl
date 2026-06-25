@testitem "fit" begin
    using FluctuationAnalysis

    @testset "exact power law is recovered" begin
        scales = [4, 8, 16, 32, 64]
        exponent = 0.75
        fluctuations = 2.0 .* float.(scales) .^ exponent
        fit = loglog_fit(scales, fluctuations)
        @test isapprox(fit.exponent, exponent; atol = 1.0e-10)
        @test isapprox(fit.intercept, log(2.0); atol = 1.0e-10)
        @test isapprox(fit.rsquared, 1.0; atol = 1.0e-12)
        @test fit.fitted_scale_range == (4, 64)
    end

    @testset "negative slope is recovered" begin
        scales = [4, 8, 16, 32, 64]
        fluctuations = 10.0 .* float.(scales) .^ (-0.5)
        @test isapprox(loglog_fit(scales, fluctuations).exponent, -0.5; atol = 1.0e-10)
    end

    @testset "noisy data gives rsquared below one" begin
        scales = [4, 8, 16, 32, 64, 128]
        fluctuations = float.(scales) .^ 0.5 .* (1 .+ 0.1 .* sin.(float.(scales)))
        @test loglog_fit(scales, fluctuations).rsquared < 1.0
    end

    @testset "fit range restricts the scales used" begin
        scales = [4, 8, 16, 32, 64]
        fluctuations = float.(scales)
        fit = loglog_fit(scales, fluctuations; fitrange = (8, 32))
        @test fit.fitted_scale_range == (8, 32)
    end

    @testset "fit range keeping exactly two scales" begin
        scales = [4, 8, 16, 32, 64]
        fit = loglog_fit(scales, float.(scales); fitrange = (16, 32))
        @test fit.fitted_scale_range == (16, 32)
    end

    @testset "coefficient of determination" begin
        observed = [1.0, 2.0, 3.0]
        @test FluctuationAnalysis.coefficient_of_determination(observed, observed) == 1.0
    end

    @testset "validation" begin
        @test_throws ArgumentError loglog_fit([4, 8], [1.0])              # length mismatch
        @test_throws ArgumentError loglog_fit([4, 8], [1.0, 2.0]; fitrange = (100, 200))
        @test_throws ArgumentError loglog_fit([4, 8], [0.0, 2.0])        # zero fluctuation
        @test_throws ArgumentError loglog_fit([4, 8], [-1.0, 2.0])       # negative
        @test_throws ArgumentError loglog_fit([4, 8], [NaN, 2.0])        # non-finite
        @test_throws ArgumentError loglog_fit([0, 8], [1.0, 2.0])        # non-positive scale
    end
end
