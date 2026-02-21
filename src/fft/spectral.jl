# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct SpectralSetup <: AbstractSetup
    sr::Int64
end

# ---------------------------------------------------------------------------- #
#                                    utils                                     #
# ---------------------------------------------------------------------------- #
# Calculate sums
sum_s = (s) -> sum(s, dims = 2)
sum_sfreq = (s, sfreq) -> sum(s .* sfreq', dims = 2)
arithmetic_mean = (s) -> sum_s(s) ./ size(s, 2)

# ---------------------------------------------------------------------------- #
#                               spectral centroid                              #
# ---------------------------------------------------------------------------- #
struct SpectralCentroid{T} <: AbstractSpectral
    spec :: Vector{<:AbstractFloat} 
    info :: SpectralSetup

    function SpectralCentroid(
        spec :: AbstractSpectrogram,    
    )
        s, sfreq = get_data(spec), get_freq(spec)
		vspec = vec(sum_sfreq(s, sfreq) ./ sum_s(s))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral crest                              #
# ---------------------------------------------------------------------------- #
struct SpectralCrest{T} <: AbstractSpectral
    spec :: Vector{<:AbstractFloat} 
    info :: SpectralSetup

    function SpectralCrest(
        spec :: AbstractSpectrogram,    
    )
        s = get_data(spec)
		peak = maximum(s, dims = 2)
		vspec = vec(peak ./ arithmetic_mean(s))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral decrease                              #
# ---------------------------------------------------------------------------- #
struct SpectralDecrease{T} <: AbstractSpectral
    spec :: Vector{<:AbstractFloat} 
    info :: SpectralSetup

    function SpectralDecrease(
        spec :: AbstractSpectrogram,    
    )
        s = get_data(spec)
		vspec = vec(sum((s[:, 2:end] .- s[:, 1]) ./ (1:size(s, 2)-1)', dims = 2) ./ sum(s[:, 2:end], dims = 2))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral entropy                               #
# ---------------------------------------------------------------------------- #
struct SpectralEntropy{T} <: AbstractSpectral
    spec :: Vector{<:AbstractFloat} 
    info :: SpectralSetup

    function SpectralEntropy(
        spec :: AbstractSpectrogram,    
    )
        s = get_data(spec)
		X = s ./ repeat(sum_s(s)', size(s, 2), 1)'
		t = replace!(-sum(X .* log2.(X), dims = 2), NaN => 0)
		vspec = vec(t ./ log2(size(s, 2)))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral flatness                              #
# ---------------------------------------------------------------------------- #
struct SpectralFlatness{T} <: AbstractSpectral
    spec :: Vector{<:AbstractFloat} 
    info :: SpectralSetup

    function SpectralFlatness(
        spec :: AbstractSpectrogram,    
    )
        s = get_data(spec)
		geometric_mean = exp.(sum(log.(s .+ eps(Float64)), dims = 2) / size(s, 2))
		vspec = vec(geometric_mean ./ arithmetic_mean(s))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                 spectral flux                                #
# ---------------------------------------------------------------------------- #
struct SpectralFlux{T} <: AbstractSpectral
    spec :: Vector{<:AbstractFloat} 
    info :: SpectralSetup

    function SpectralFlux(
        spec :: AbstractSpectrogram;
		p :: Int64=2
    )
        s = get_data(spec)
		initial_condition = s[1, :]
		temp = diff(hcat(initial_condition, s'), dims = 2)
		vspec = [LinearAlgebra.norm(temp[:, i], p) for i in axes(temp, 2)]
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    get_data(s::Spectral) -> Vector

Get the Spectral data vector.
"""
@inline get_data(s::AbstractSpectral)  = s.spec

"""
    get_setup(s::Spectral) -> MelSpecSetup

Get the configuration metadata for the Spectral spectrogram.
"""
@inline get_setup(s::AbstractSpectral) = s.info











# function get_spectrum(setup::AudioSetup, data::AudioData)
# 	setup.spectral_spectrum == :mel && return data.mel_spectrogram', data.mel_frequencies
# 	setup.spectral_spectrum == :lin && return data.lin_spectrogram', data.lin_frequencies
# 	error("Unknown spectral spectrum")
# end

# function spectral_kurtosis(
# 	s::AbstractArray{Float64},
# 	data::AudioData,
# 	# sum_x1::Vector{Float64},
# 	# centroid::Vector{Float64},
# 	higher_moment_tmp::AbstractArray{Float64},
# 	higher_moment_denom::Vector{Float64},
# 	higher_momement_num::AbstractArray{Float64},
# )
# 	# calculate centroid
# 	# centroid = vec(sum(s .* setup.fft_frequencies, dims=1) ./ sum(s, dims=1))
# 	# # calculate spread
# 	# temp = (setup.fft_frequencies .- centroid') .^ 2
# 	# spread = sqrt.(sum((temp) .* s, dims=1) ./ sum(s, dims=1))
# 	# # calculate kurtosis
# 	# data.spectral_kurtosis = vec(real(sum(((temp .^ 2) .* s), dims=1) ./ ((spread .^ 4) .* sum(s, dims=1))))
# 	data.spectral_kurtosis = vec(sum(higher_momement_num .* higher_moment_tmp, dims = 1) ./ (higher_moment_denom .* data.spectral_spread)')
# end

# function spectral_rolloff(s::AbstractArray{Float64}, data::AudioData, fft_frequencies::Vector{Float64})
# 	# calculate rolloff point
# 	threshold = 0.95
# 	c = cumsum(s, dims = 1)
# 	d = c[end, :] * threshold
# 	idx = zeros(size(d))
# 	for i in eachindex(d)
# 		idx[i] = findfirst(real(c[:, i]) .>= real(d[i]))
# 	end
# 	data.spectral_rolloff = fft_frequencies[Int.(idx)]
# end

# function spectral_skewness(
# 	s::AbstractArray{Float64},
# 	data::AudioData,
# 	sum_x1::Vector{Float64},
# 	centroid::Vector{Float64},
# 	higher_moment_denom::Vector{Float64},
# 	higher_momement_num::AbstractArray{Float64},
# )
# 	# # calculate centroid
# 	# centroid = vec(sum(s .* setup.fft_frequencies, dims=1) ./ sum(s, dims=1))
# 	# # calculate spread
# 	# temp = (setup.fft_frequencies .- centroid')
# 	# spread = sqrt.(sum((temp .^ 2) .* s, dims=1) ./ sum(s, dims=1))
# 	# # calculate skewness
# 	# data.spectral_skewness = vec(real(sum((temp .^ 3) .* s, dims=1) ./ ((spread .^ 3) .* sum(s, dims=1))))
# 	data.spectral_skewness = vec(sum(higher_momement_num, dims = 1) ./ higher_moment_denom')
# end

# function spectral_slope(
# 	s::AbstractArray{Float64},
# 	data::AudioData,
# 	sum_x1::Vector{Float64},
# 	arithmetic_mean::Vector{Float64},
# 	fft_frequencies::Vector{Float64},
# )
# 	# calculate slope
# 	f_minus_mu_f = fft_frequencies .- sum(fft_frequencies, dims = 1) ./ size(s, 1)
# 	X_minus_mu_X = s .- arithmetic_mean'
# 	data.spectral_slope = vec(real(sum(X_minus_mu_X .* f_minus_mu_f, dims = 1) ./ sum(f_minus_mu_f .^ 2)))
# end

# function _get_spread(x::AbstractArray{Float64}, setup::AudioSetup)

# end

# ################################################################################
# #                                    main                                      #
# ################################################################################
# function get_spectrals!(
# 	setup::AudioSetup,
# 	data::AudioData,
# )
# 	s, freq = get_spectrum(setup, data)

# 	# common data
# 	size_x1 = size(s, 1)
# 	sum_x1 = vec(sum(s, dims = 1))
# 	arithmetic_mean = sum_x1 ./ size_x1
# 	data.spectral_centroid = vec(sum(s .* freq, dims = 1) ./ sum_x1')
# 	data.spectral_centroid = replace!(data.spectral_centroid, NaN => 0)
# 	higher_moment_tmp = freq .- data.spectral_centroid'
# 	data.spectral_spread = vec(sqrt.(sum((higher_moment_tmp .^ 2) .* s, dims = 1) ./ sum_x1'))
# 	higher_moment_denom = (data.spectral_spread .^ 3) .* sum_x1
# 	higher_momement_num = (higher_moment_tmp .^ 3) .* s

# 	spectral_crest(s, data, sum_x1, arithmetic_mean)
# 	spectral_entropy(s, data, sum_x1)
# 	spectral_flatness(s, data, sum_x1, arithmetic_mean)
# 	spectral_flux(s, data)
# 	spectral_kurtosis(s, data, higher_moment_tmp, higher_moment_denom, higher_momement_num)
# 	spectral_rolloff(s, data, freq)
# 	spectral_skewness(s, data, sum_x1, data.spectral_centroid, higher_moment_denom, higher_momement_num)
# 	spectral_decrease(s, data)
# 	spectral_slope(s, data, sum_x1, arithmetic_mean, freq)
# end

# # ---------------------------------------------------------------------------------- #
# #                                     utilities                                      #
# # ---------------------------------------------------------------------------------- #
# # Calculate sums
# sum_s = (s) -> sum(s, dims = 1)
# sum_sfreq = (s, sfreq) -> sum(s .* sfreq, dims = 1)


# # ---------------------------------------------------------------------------------- #
# #                                  spectral spread                                   #
# # ---------------------------------------------------------------------------------- #
# function _get_spec_spread(
# 	s::AbstractArray{Float64},
# 	sfreq::AbstractVector{Float64},
# )
# 	centroid = _get_spec_centroid(s, sfreq)
# 	vec(sqrt.(sum_s(s .* (sfreq .- centroid').^2) ./ sum_s(s)))
# end