"""
Spectral Centroid:
The spectral centroid represents the "center of gravity" of the spectrum. 
It is used as an indication of brightness and is commonly used in music analysis and genre classification. 
For example, observe the jumps in the centroid corresponding to high hat hits in the audio file.

Spectral Entropy:
Spectral entropy quantifies the uncertainty or randomness in the spectral distribution.
It has been used successfully in voiced/unvoiced decisions for automatic speech recognition. 
Because entropy is a measure of disorder, regions of voiced speech have lower entropy compared to regions of unvoiced speech.

Spectral Flatness:
Spectral flatness measures the balance between tonal and noisy components in the speech signal. 
It indicates the presence of harmonic or noisy characteristics.

Spectral Kurtosis:
Spectral kurtosis indicates the “tailedness” of the spectrum. 
It offers insights into the presence of abrupt spectral changes and harmonic components in speech.

Spectral Roll-Off:
Spectral roll-off frequency is the point below which a certain percentage of the total spectral energy lies. 
It helps distinguish between voiced and unvoiced speech sounds.

Spectral Slope:
Spectral slope measures the rate of change in the spectral distribution.
It has been used extensively in speech analysis, particularly in modeling speaker stress.
The slope is directly related to the resonant characteristics of the vocal folds and has also been applied to speaker identification. 
it is a socially important aspect of timbre. 
Spectral slope discrimination has been shown to occur in early childhood development. 
It is most pronounced when the energy in the lower formants is much greater than the energy in the higher formants.

Spectral Spread:
Spectral spread characterizes the width or bandwidth of the frequency distribution in the spectrum. 
It indicates the extent of frequency variations present in the speech signal.

Spectral Skewness:
Spectral skewness measures the asymmetry of the frequency distribution. 
It reveals whether the speech sound is dominated by higher or lower frequencies.
"""
# ---------------------------------------------------------------------------- #
#                             spectral features                                #
# ---------------------------------------------------------------------------- #
struct SpectralSetup
    freqrange::Tuple{Int, Int}
end

struct SpectralData{T<:AbstractFloat}
	centroid::AbstractVector{T}
	crest::AbstractVector{T}
	decrease::AbstractVector{T}
	entropy::AbstractVector{T}
	flatness::AbstractVector{T}
	flux::AbstractVector{T}
	kurtosis::AbstractVector{T}
	rolloff::AbstractVector{T}
	skewness::AbstractVector{T}
	slope::AbstractVector{T}
	spread::AbstractVector{T}
end

struct Spectral
    sr::Int
    setup::SpectralSetup
    data::SpectralData
end

# ---------------------------------------------------------------------------------- #
#                                     utilities                                      #
# ---------------------------------------------------------------------------------- #
# Calculate sums
sum_s = (s) -> sum(s, dims = 1)
sum_sfreq = (s, sfreq) -> sum(s .* sfreq, dims = 1)

# ---------------------------------------------------------------------------------- #
#                             single spectral functions                              #
# ---------------------------------------------------------------------------------- #
function spectral_crest(s::AbstractArray{T}, arithmetic_mean::Vector{T}) where T <: AbstractFloat
	# calculate spectral peak
	peak = maximum(s, dims = 1)
	# calculate spectral crest
	return vec(peak ./ arithmetic_mean')
end

function spectral_decrease(s::AbstractArray{<:AbstractFloat})
	# calculate decrease
	return vec(real(sum((s[2:end, :] .- s[1, :]') ./ (1:size(s, 1)-1), dims = 1) ./ sum(s[2:end, :], dims = 1)))
end

function spectral_entropy(s::AbstractArray{T}, sum_x1::Vector{T}) where T <: AbstractFloat
	# calculate entropy
	X = s ./ repeat(sum_x1', size(s, 1), 1)
	t = replace!(-sum(X .* log2.(X), dims = 1), NaN => 0)
	return vec(t ./ log2(size(s, 1)))
end

function spectral_flatness(s::AbstractArray{T}, arithmetic_mean::Vector{T}) where T <: AbstractFloat
	# calculate geometric mean, arithmetic mean, and flatness
	geometric_mean = exp.(sum(log.(s .+ eps(T)), dims = 1) / size(s, 1))
	return vec(geometric_mean ./ arithmetic_mean')
end

function spectral_flux(s::AbstractArray{<:AbstractFloat})
	initial_condition = s[:, 1]
	# calculate flux
	temp = diff(hcat(initial_condition, s), dims = 2)
	fl = Float64[]
	for i in axes(temp, 2)
		append!(fl, norm(temp[:, i]))
	end
	return fl
end

function spectral_kurtosis(spread::AbstractVector{T}, higher_moment_tmp::AbstractArray{T}, higher_moment_denom::Vector{T}, higher_momement_num::AbstractArray{T}) where T <: AbstractFloat
	return vec(sum(higher_momement_num .* higher_moment_tmp, dims = 1) ./ (higher_moment_denom .* spread)')
end

function spectral_rolloff(s::AbstractArray{T}, fft_frequencies::AbstractVector{T}) where T <: AbstractFloat
	# calculate rolloff point
	threshold = 0.95
	c = cumsum(s, dims = 1)
	d = c[end, :] * threshold
	idx = zeros(size(d))
	for i in eachindex(d)
		idx[i] = findfirst(real(c[:, i]) .>= real(d[i]))
	end
	return fft_frequencies[Int.(idx)]
end

function spectral_skewness(higher_moment_denom::Vector{T}, higher_momement_num::AbstractArray{T}) where T <: AbstractFloat
	return vec(sum(higher_momement_num, dims = 1) ./ higher_moment_denom')
end

function spectral_slope(s::AbstractArray{T}, arithmetic_mean::AbstractVector{T}, fft_frequencies::AbstractVector{T}) where T <: AbstractFloat
	# calculate slope
	f_minus_mu_f = fft_frequencies .- sum(fft_frequencies, dims = 1) ./ size(s, 1)
	X_minus_mu_X = s .- arithmetic_mean'
	return vec(real(sum(X_minus_mu_X .* f_minus_mu_f, dims = 1) ./ sum(f_minus_mu_f .^ 2)))
end

# ---------------------------------------------------------------------------- #
#                             spectral features                                #
# ---------------------------------------------------------------------------- #
function _get_spectrals(
	s::AbstractArray{T},
	sr::Int64;
	freq::AbstractVector{T},
	freqrange = (0, round(Int, sr / 2))
) where {T <: AbstractFloat}
	# trim to desired frequency range
	s_range = findall(freqrange[1] .<= freq .<= freqrange[2])
	s, freq = s[s_range, :], freq[s_range]

	# common datas
	size_x1 = size(s, 1)
	sum_x1 = vec(sum(s, dims=1))
	arithmetic_mean = sum_x1 ./ size_x1
	centroid = replace!(vec(sum(s .* freq, dims = 1) ./ sum_x1'), NaN => 0)
	higher_moment_tmp = freq .- centroid'
	spread = vec(sqrt.(sum((higher_moment_tmp .^ 2) .* s, dims = 1) ./ sum_x1'))
	higher_moment_denom = (spread .^ 3) .* sum_x1
	higher_momement_num = (higher_moment_tmp .^ 3) .* s

	crest = spectral_crest(s, arithmetic_mean)
	entropy = spectral_entropy(s, sum_x1)
	flatness = spectral_flatness(s, arithmetic_mean)
	flux = spectral_flux(s)
	kurtosis = spectral_kurtosis(spread, higher_moment_tmp, higher_moment_denom, higher_momement_num)
	rolloff = spectral_rolloff(s, freq)
	skewness = spectral_skewness(higher_moment_denom, higher_momement_num)
	decrease = spectral_decrease(s)
	slope = spectral_slope(s, arithmetic_mean, freq)

	Spectral(sr, SpectralSetup(freqrange), SpectralData(centroid, crest, decrease, entropy, flatness, flux, kurtosis, rolloff, skewness, slope, spread))
end

function Base.show(io::IO, s::Spectral)
    println(io, "Spectral Features:")
    println(io, "  Sample Rate: $(s.sr) Hz")
    println(io, "  Frequency Range: $(s.setup.freqrange) Hz")
    
    features = [
        ("Centroid", s.data.centroid),
        ("Crest", s.data.crest),
        ("Decrease", s.data.decrease),
        ("Entropy", s.data.entropy),
        ("Flatness", s.data.flatness),
        ("Flux", s.data.flux),
        ("Kurtosis", s.data.kurtosis),
        ("Rolloff", s.data.rolloff),
        ("Skewness", s.data.skewness),
        ("Slope", s.data.slope),
        ("Spread", s.data.spread)
    ]
    
    for (name, feature) in features
        status = isnothing(feature) ? "Not computed" : "Computed"
        println(io, "  $name: $status")
    end
end

using Plots

function Base.display(s::Spectral)
    features = [
        ("Centroid", s.data.centroid),
        ("Crest", s.data.crest),
        ("Decrease", s.data.decrease),
        ("Entropy", s.data.entropy),
        ("Flatness", s.data.flatness),
        ("Flux", s.data.flux),
        ("Kurtosis", s.data.kurtosis),
        ("Rolloff", s.data.rolloff),
        ("Skewness", s.data.skewness),
        ("Slope", s.data.slope),
        ("Spread", s.data.spread)
    ]

    n_features = count(!isnothing, last.(features))
    plots = []

    for (name, feature) in features
        if !isnothing(feature)
            push!(plots, plot(feature, title=name, label=""))
        end
    end

    plot(plots..., layout=(n_features, 1), size=(800, 300 * n_features), link=:x)
end


function get_spectrals(source::Union{Stft, Cwt}; kwargs...)
	_get_spectrals(source.data.spec, source.sr; freq=source.data.freq, kwargs...)
end
