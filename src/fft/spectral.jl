# using LinearAlgebra

# include("../signalDataStructure.jl")
# include("fft.jl")

function get_spectrum(data::signal_data, setup::signal_setup)
    setup.spectral_spectrum == :mel && return data.mel_spectrogram', setup.mel_frequencies
    setup.spectral_spectrum == :linear && return data.lin_spectrogram', setup.lin_frequencies
    error("Unknown spectral spectrum")
end

function spectral_crest(
    s::AbstractArray{Float64},
    data::signal_data,
    sum_x1::Vector{Float64},
    arithmetic_mean::Vector{Float64}
)
    # # calculate spectral mean
    # m = sum(real(s), dims=1) ./ size(s, 1)
    # calculate spectral peak
    peak = maximum(s, dims=1)
    # calculate spectral crest
    data.spectral_crest = vec(peak ./ arithmetic_mean')
end

function spectral_decrease(s::AbstractArray{Float64}, data::signal_data)
    # calculate decrease
    data.spectral_decrease = vec(real(sum((s[2:end, :] .- s[1, :]') ./ (1:size(s, 1)-1), dims=1) ./ sum(s[2:end, :], dims=1)))
end

function spectral_entropy(
    s::AbstractArray{Float64},
    data::signal_data,
    sum_x1::Vector{Float64}
)
    # calculate entropy
    X = s ./ repeat(sum_x1', size(s, 1), 1)
    t = replace!(-sum(X .* log2.(X), dims=1), NaN => 0)
    data.spectral_entropy = vec(t ./ log2(size(s, 1)))
end

function spectral_flatness(
    s::AbstractArray{Float64},
    data::signal_data,
    sum_x1::Vector{Float64},
    arithmetic_mean::Vector{Float64}
)
    # calculate geometric mean, arithmetic mean, and flatness
    geometric_mean = exp.(sum(log.(s .+ eps(Float64)), dims=1) / size(s, 1))
    data.spectral_flatness = vec(geometric_mean ./ arithmetic_mean')
end

function spectral_flux(s::AbstractArray{Float64}, data::signal_data)
    initial_condition = s[:, 1]
    # calculate flux
    temp = diff(hcat(initial_condition, s), dims=2)
    fl = []
    for i in axes(temp, 2)
        append!(fl, norm(temp[:, i]))
    end
    data.spectral_flux = fl
end

function spectral_kurtosis(
    s::AbstractArray{Float64},
    data::signal_data,
    # sum_x1::Vector{Float64},
    # centroid::Vector{Float64},
    higher_moment_tmp::AbstractArray{Float64},
    higher_moment_denom::Vector{Float64},
    higher_momement_num::AbstractArray{Float64}
)
    # calculate centroid
    # centroid = vec(sum(s .* setup.fft_frequencies, dims=1) ./ sum(s, dims=1))
    # # calculate spread
    # temp = (setup.fft_frequencies .- centroid') .^ 2
    # spread = sqrt.(sum((temp) .* s, dims=1) ./ sum(s, dims=1))
    # # calculate kurtosis
    # data.spectral_kurtosis = vec(real(sum(((temp .^ 2) .* s), dims=1) ./ ((spread .^ 4) .* sum(s, dims=1))))
    data.spectral_kurtosis = vec(sum(higher_momement_num .* higher_moment_tmp, dims=1) ./ (higher_moment_denom .* data.spectral_spread)')
end

function spectral_rolloff(s::AbstractArray{Float64}, data::signal_data, fft_frequencies::Vector{Float64})
    # calculate rolloff point
    threshold = 0.95
    c = cumsum(s, dims=1)
    d = c[end, :] * threshold
    idx = zeros(size(d))
    for i in eachindex(d)
        idx[i] = findfirst(real(c[:, i]) .>= real(d[i]))
    end
    data.spectral_rolloff = fft_frequencies[Int.(idx)]
end

function spectral_skewness(
    s::AbstractArray{Float64},
    data::signal_data,
    sum_x1::Vector{Float64},
    centroid::Vector{Float64},
    higher_moment_denom::Vector{Float64},
    higher_momement_num::AbstractArray{Float64}
)
    # # calculate centroid
    # centroid = vec(sum(s .* setup.fft_frequencies, dims=1) ./ sum(s, dims=1))
    # # calculate spread
    # temp = (setup.fft_frequencies .- centroid')
    # spread = sqrt.(sum((temp .^ 2) .* s, dims=1) ./ sum(s, dims=1))
    # # calculate skewness
    # data.spectral_skewness = vec(real(sum((temp .^ 3) .* s, dims=1) ./ ((spread .^ 3) .* sum(s, dims=1))))
    data.spectral_skewness = vec(sum(higher_momement_num, dims=1) ./ higher_moment_denom')
end

function spectral_slope(
    s::AbstractArray{Float64},
    data::signal_data,
    sum_x1::Vector{Float64},
    arithmetic_mean::Vector{Float64},
    fft_frequencies::Vector{Float64}
)
    # calculate slope
    f_minus_mu_f = fft_frequencies .- sum(fft_frequencies, dims=1) ./ size(s, 1)
    X_minus_mu_X = s .- arithmetic_mean'
    data.spectral_slope = vec(real(sum(X_minus_mu_X .* f_minus_mu_f, dims=1) ./ sum(f_minus_mu_f .^ 2)))
end

# function spectral_features(
#     x::AbstractArray{T},
#     sr::Int64;
#     window_length::Int64=Int(round(0.03 * sr)),
#     overlap_length::Int64=Int(round(0.02 * sr)),
#     frequency_range::Vector{Int64}=[0, Int(round(sr / 2))],
#     window_type::Symbol=:hann,
#     # windowNormalization::Bool=true,
#     # oneSided::Bool=true
# ) where {T<:AbstractFloat}

#     # options and data structures definition
#     options = signal_setup(
#         sr=sr,
#         window_length=window_length,
#         overlap_length=overlap_length,
#         frequency_range=Float64.(frequency_range),
#         window_type=window_type,
#         spectrum_type=:power,
#         # windowNormalization=windowNormalization,
#         # oneSided=oneSided
#     )

#     data = signal_data(
#         x=Float64.(x)
#     )

#     takeFFT(data, options)
#     spectral_features(data, options)
# end

function spectral_features(data::signal_data, setup::signal_setup)
    s, freq = get_spectrum(data, setup)

    # common data
    size_x1 = size(s, 1)
    sum_x1 = vec(sum(s, dims=1))
    arithmetic_mean = sum_x1 ./ size_x1
    data.spectral_centroid = vec(sum(s .* freq, dims=1) ./ sum_x1')
    data.spectral_centroid = replace!(data.spectral_centroid, NaN => 0)
    higher_moment_tmp = freq .- data.spectral_centroid'
    data.spectral_spread = vec(sqrt.(sum((higher_moment_tmp .^ 2) .* s, dims=1) ./ sum_x1'))
    higher_moment_denom = (data.spectral_spread .^ 3) .* sum_x1
    higher_momement_num = (higher_moment_tmp .^ 3) .* s

    spectral_crest(s, data, sum_x1, arithmetic_mean)
    spectral_entropy(s, data, sum_x1)
    spectral_flatness(s, data, sum_x1, arithmetic_mean)
    spectral_flux(s, data)
    spectral_kurtosis(s, data, higher_moment_tmp, higher_moment_denom, higher_momement_num)
    spectral_rolloff(s, data, freq)
    spectral_skewness(s, data, sum_x1, data.spectral_centroid, higher_moment_denom, higher_momement_num)
    spectral_decrease(s, data)
    spectral_slope(s, data, sum_x1, arithmetic_mean, freq)
end

function spectral_spread(
    x::AbstractArray{T},
    sr::Int64;
    fft_length::Int64=256,
    window_type::Vector{Symbol}=[:hann, :periodic],
    window_length::Int64=Int(round(0.03 * sr)),
    overlap_length::Int64=Int(round(0.02 * sr)),
    window_norm::Bool=true,
    frequency_range::Vector{Int64}=[0, Int(floor(sr / 2))],
    spectrum_type::Symbol=:magnitude
) where {T<:AbstractFloat}
    # setup and data structures definition    
    setup_spread = signal_setup(
        sr=sr,
        fft_length=fft_length,
        window_type=window_type,
        window_length=window_length,
        overlap_length=overlap_length,
        window_norm=window_norm,
        frequency_range=frequency_range,
        spectrum_type=spectrum_type
    )

    data_spread = signal_data(
        x=Float64.(x)
    )

    takeFFT(data_spread, setup_spread)
    lin_spectrogram(data_spread, setup_spread)

    s, freq = data_spread.lin_spectrogram', setup_spread.lin_frequencies

    sum_x1 = vec(sum(s, dims=1))
    data_spread.spectral_centroid = vec(sum(s .* freq, dims=1) ./ sum_x1')
    data_spread.spectral_centroid = replace!(data_spread.spectral_centroid, NaN => 0)
    higher_moment_tmp = freq .- data_spread.spectral_centroid'

    data_spread.spectral_spread = vec(sqrt.(sum((higher_moment_tmp .^ 2) .* s, dims=1) ./ sum_x1'))
    
    return data_spread.spectral_spread
end

