@testitem "segmentation" begin
    using FluctuationAnalysis
    segment_views = FluctuationAnalysis.segment_views
    profile = collect(1.0:10.0)

    @testset "non-overlapping forward segmentation" begin
        segments = segment_views(profile, 5; overlap = false, bidirectional = false)
        @test length(segments) == 2
        @test collect(segments[1]) == [1.0, 2.0, 3.0, 4.0, 5.0]
        @test collect(segments[2]) == [6.0, 7.0, 8.0, 9.0, 10.0]
    end

    @testset "bidirectional covers the remainder" begin
        segments = segment_views(profile, 4; overlap = false, bidirectional = true)
        @test length(segments) == 4
        @test collect(segments[end]) == [7.0, 8.0, 9.0, 10.0]
    end

    @testset "bidirectional doubles the count when the scale divides evenly" begin
        segments = segment_views(profile, 5; overlap = false, bidirectional = true)
        @test length(segments) == 4   # two forward + two from the end (same coverage)
    end

    @testset "overlapping sliding window" begin
        segments = segment_views(profile, 4; overlap = true)
        @test length(segments) == length(profile) - 4 + 1
        @test collect(segments[1]) == [1.0, 2.0, 3.0, 4.0]
        @test collect(segments[end]) == [7.0, 8.0, 9.0, 10.0]
    end

    @testset "smallest scale" begin
        segments = segment_views(profile, 2; overlap = false, bidirectional = false)
        @test length(segments) == 5
        @test collect(segments[1]) == [1.0, 2.0]
    end

    @testset "scale equal to the profile length yields one segment" begin
        segments = segment_views(profile, 10; overlap = false, bidirectional = false)
        @test length(segments) == 1
        @test collect(segments[1]) == profile
        @test length(segment_views(profile, 10; overlap = true)) == 1
    end

    @testset "segments are views into the original profile" begin
        segments = segment_views(profile, 5; overlap = false, bidirectional = false)
        @test segments[1] isa SubArray
    end

    @testset "validation" begin
        @test_throws ArgumentError segment_views(profile, 1)
        @test_throws ArgumentError segment_views(profile, 11)
    end
end
