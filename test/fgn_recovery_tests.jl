# Scientific recovery tests: the estimators must recover an *intermediate* Hurst
# exponent (not just the H = 0.5 / H = 1.5 special points) from a signal with known
# long-range correlations. The signal is exact fractional Gaussian noise generated
# by the Hosking (Durbin–Levinson) recursion — the same generator validated to
# machine precision in paper/src/ReferenceData.jl, ported here so the package
# tests stay dependency-free.

@testsetup module FractionalNoise
using Random: AbstractRNG, MersenneTwister, randn
export fractional_gaussian_noise

fgn_autocovariance(lag::Integer, hurst::Real) =
    0.5 * (abs(lag + 1)^(2hurst) - 2 * abs(lag)^(2hurst) + abs(lag - 1)^(2hurst))

"Exact fGn of Hurst `hurst` (unit variance) via the Hosking recursion."
function fractional_gaussian_noise(n_samples::Integer, hurst::Real; rng::AbstractRNG)
    output = Vector{Float64}(undef, n_samples)
    phi = zeros(Float64, n_samples)
    previous = zeros(Float64, n_samples)
    output[1] = randn(rng)
    variance = 1.0
    n_samples == 1 && return output
    phi[1] = fgn_autocovariance(1, hurst)
    variance *= 1 - phi[1]^2
    output[2] = phi[1] * output[1] + sqrt(variance) * randn(rng)
    for step in 2:(n_samples - 1)
        copyto!(previous, 1, phi, 1, step - 1)
        numerator = fgn_autocovariance(step, hurst)
        for index in 1:(step - 1)
            numerator -= previous[index] * fgn_autocovariance(step - index, hurst)
        end
        reflection = numerator / variance
        for index in 1:(step - 1)
            phi[index] = previous[index] - reflection * previous[step - index]
        end
        phi[step] = reflection
        variance *= 1 - reflection^2
        mean_estimate = 0.0
        for index in 1:step
            mean_estimate += phi[index] * output[step + 1 - index]
        end
        output[step + 1] = mean_estimate + sqrt(variance) * randn(rng)
    end
    return output
end
end

@testitem "dfa: intermediate Hurst recovery from fGn" setup = [FractionalNoise] begin
    using FluctuationAnalysis
    using Random: MersenneTwister
    using Statistics: mean

    rng = MersenneTwister(11)
    scales = logarithmic_scales(
        4096; minimum_scale = 16, maximum_scale = 1024, scale_count = 20
    )
    for hurst_value in (0.3, 0.6, 0.8)
        estimates = [
            scaling_exponent(
                    dfa(
                        FractionalNoise.fractional_gaussian_noise(
                            4096, hurst_value; rng = rng
                        );
                        scales = scales,
                    ),
                )
                for _ in 1:4
        ]
        @test isapprox(mean(estimates), hurst_value; atol = 0.06)
    end
end

@testitem "hurst: intermediate Hurst recovery from fGn" setup = [FractionalNoise] begin
    using FluctuationAnalysis
    using Random: MersenneTwister
    using Statistics: mean

    rng = MersenneTwister(13)
    scales = logarithmic_scales(
        4096; minimum_scale = 16, maximum_scale = 1024, scale_count = 20
    )
    for hurst_value in (0.3, 0.7)
        series = [
            FractionalNoise.fractional_gaussian_noise(4096, hurst_value; rng = rng)
                for _ in 1:4
        ]
        dfa_hurst = mean(
            hurst_exponent(
                    hurst(value, DetrendedFluctuationHurst(); scales = scales)
                ) for value in series
        )
        @test isapprox(dfa_hurst, hurst_value; atol = 0.06)
        # Classic R/S carries a known finite-sample bias, so it gets a looser bound.
        rs_hurst = mean(
            hurst_exponent(hurst(value, RescaledRangeHurst(); scales = scales))
                for value in series
        )
        @test isapprox(rs_hurst, hurst_value; atol = 0.13)
    end
end

@testitem "dcca: coupled fGn cross-exponent" setup = [FractionalNoise] begin
    using FluctuationAnalysis
    using Random: MersenneTwister
    using Statistics: mean

    rng = MersenneTwister(17)
    hurst_value = 0.7
    scales = logarithmic_scales(
        4096; minimum_scale = 16, maximum_scale = 1024, scale_count = 20
    )
    lambdas = Float64[]
    for _ in 1:4
        shared = FractionalNoise.fractional_gaussian_noise(4096, hurst_value; rng = rng)
        first_noise = FractionalNoise.fractional_gaussian_noise(
            4096, hurst_value; rng = rng
        )
        second_noise = FractionalNoise.fractional_gaussian_noise(
            4096, hurst_value; rng = rng
        )
        first_series = sqrt(0.5) .* shared .+ sqrt(0.5) .* first_noise
        second_series = sqrt(0.5) .* shared .+ sqrt(0.5) .* second_noise
        push!(lambdas, scaling_exponent(dcca(first_series, second_series; scales = scales)))
    end
    # Every component shares Hurst H, so the cross-correlation exponent λ ≈ H.
    @test isapprox(mean(lambdas), hurst_value; atol = 0.07)
end
