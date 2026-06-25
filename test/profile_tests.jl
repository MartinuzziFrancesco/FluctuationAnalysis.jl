@testitem "profile" begin
    using FluctuationAnalysis
    using Statistics

    @testset "demeaned profile ends at zero" begin
        series = [1.0, 2.0, 3.0, 4.0]
        result = integrated_profile(series)
        @test length(result) == length(series)
        @test isapprox(result[end], 0.0; atol = 1.0e-12)
    end

    @testset "known cumulative sum" begin
        series = [2.0, 4.0, 6.0]
        @test integrated_profile(series; demean = false) == [2.0, 6.0, 12.0]
    end

    @testset "demean removes exactly the mean" begin
        series = [5.0, 1.0, 4.0, 2.0]
        @test integrated_profile(series) == cumsum(series .- mean(series))
    end

    @testset "integer input is promoted to Float64" begin
        result = integrated_profile([1, 2, 3]; demean = false)
        @test eltype(result) == Float64
        @test result == [1.0, 3.0, 6.0]
    end

    @testset "single-element profile" begin
        @test integrated_profile([7.0]) == [0.0]
        @test integrated_profile([7.0]; demean = false) == [7.0]
    end

    @testset "result is always a one-based Vector" begin
        result = integrated_profile([3.0, -1.0, 2.0, 6.0])
        @test result isa Vector{Float64}
        @test firstindex(result) == 1
    end

    @testset "empty input is rejected" begin
        @test_throws ArgumentError integrated_profile(Float64[])
    end
end
