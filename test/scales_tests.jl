@testitem "scales" begin
    using FluctuationAnalysis

    @testset "scales are sorted, distinct, and bounded" begin
        scales = logarithmic_scales(1000; minimum_scale = 4, maximum_scale = 200, scale_count = 20)
        @test issorted(scales)
        @test allunique(scales)
        @test minimum(scales) >= 4
        @test maximum(scales) <= 200
        @test eltype(scales) == Int
    end

    @testset "rounding collisions are deduplicated" begin
        # A wide count over a narrow range must collapse to fewer distinct scales.
        scales = logarithmic_scales(1000; minimum_scale = 4, maximum_scale = 8, scale_count = 20)
        @test allunique(scales)
        @test length(scales) <= 5
    end

    @testset "endpoints are honoured" begin
        scales = logarithmic_scales(4096; minimum_scale = 8, maximum_scale = 512, scale_count = 12)
        @test first(scales) == 8
        @test last(scales) == 512
    end

    @testset "base changes the spacing" begin
        coarse = logarithmic_scales(10_000; minimum_scale = 4, maximum_scale = 2000, base = 2.0)
        fine = logarithmic_scales(10_000; minimum_scale = 4, maximum_scale = 2000, base = 1.2)
        @test length(fine) >= length(coarse)
    end

    @testset "validation" begin
        @test_throws ArgumentError logarithmic_scales(5)
        @test_throws ArgumentError logarithmic_scales(1000; minimum_scale = 1)
        @test_throws ArgumentError logarithmic_scales(1000; minimum_scale = 50, maximum_scale = 40)
        @test_throws ArgumentError logarithmic_scales(1000; maximum_scale = 2000)
        @test_throws ArgumentError logarithmic_scales(1000; scale_count = 1)
        @test_throws ArgumentError logarithmic_scales(1000; base = 1.0)
    end
end
