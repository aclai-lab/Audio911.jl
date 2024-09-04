# ---------------------------------------------------------------------------- #
#                             spectral features                                #
# ---------------------------------------------------------------------------- #
struct SpectralSetup
    freqrange::Tuple{Int64, Int64}
end

struct SpectralData
	centroid::AbstractVector{Float64}
	crest::AbstractVector{Float64}
	decrease::AbstractVector{Float64}
	entropy::AbstractVector{Float64}
	flatness::AbstractVector{Float64}
	flux::AbstractVector{Float64}
	kurtosis::AbstractVector{Float64}
	rolloff::AbstractVector{Float64}
	skewness::AbstractVector{Float64}
	slope::AbstractVector{Float64}
	spread::AbstractVector{Float64}
end

struct Spectral
    sr::Int64
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
function spectral_crest(s::AbstractArray{Float64}, arithmetic_mean::Vector{Float64})
	# calculate spectral peak
	peak = maximum(s, dims = 1)
	# calculate spectral crest
	return vec(peak ./ arithmetic_mean')
end

function spectral_decrease(s::AbstractArray{Float64})
	# calculate decrease
	return vec(real(sum((s[2:end, :] .- s[1, :]') ./ (1:size(s, 1)-1), dims = 1) ./ sum(s[2:end, :], dims = 1)))
end

function spectral_entropy(s::AbstractArray{Float64}, sum_x1::Vector{Float64})
	# calculate entropy
	X = s ./ repeat(sum_x1', size(s, 1), 1)
	t = replace!(-sum(X .* log2.(X), dims = 1), NaN => 0)
	return vec(t ./ log2(size(s, 1)))
end

function spectral_flatness(s::AbstractArray{Float64}, arithmetic_mean::Vector{Float64})
	# calculate geometric mean, arithmetic mean, and flatness
	geometric_mean = exp.(sum(log.(s .+ eps(Float64)), dims = 1) / size(s, 1))
	return vec(geometric_mean ./ arithmetic_mean')
end

function spectral_flux(s::AbstractArray{Float64})
	initial_condition = s[:, 1]
	# calculate flux
	temp = diff(hcat(initial_condition, s), dims = 2)
	fl = []
	for i in axes(temp, 2)
		append!(fl, norm(temp[:, i]))
	end
	return fl
end

function spectral_kurtosis(spread::AbstractVector{Float64}, higher_moment_tmp::AbstractArray{Float64}, higher_moment_denom::Vector{Float64}, higher_momement_num::AbstractArray{Float64})
	return vec(sum(higher_momement_num .* higher_moment_tmp, dims = 1) ./ (higher_moment_denom .* spread)')
end

function spectral_rolloff(s::AbstractArray{Float64}, fft_frequencies::AbstractVector{Float64})
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

function spectral_skewness(higher_moment_denom::Vector{Float64}, higher_momement_num::AbstractArray{Float64})
	return vec(sum(higher_momement_num, dims = 1) ./ higher_moment_denom')
end

function spectral_slope(s::AbstractArray{Float64}, arithmetic_mean::AbstractVector{Float64}, fft_frequencies::AbstractVector{Float64})
	# calculate slope
	f_minus_mu_f = fft_frequencies .- sum(fft_frequencies, dims = 1) ./ size(s, 1)
	X_minus_mu_X = s .- arithmetic_mean'
	return vec(real(sum(X_minus_mu_X .* f_minus_mu_f, dims = 1) ./ sum(f_minus_mu_f .^ 2)))
end

# ---------------------------------------------------------------------------- #
#                             spectral features                                #
# ---------------------------------------------------------------------------- #
function _get_spectrals(;
	s::AbstractArray{Float64},
	sr::Int64,
	freq::AbstractVector{Float64},
	freqrange = (0, round(Int, sr / 2))
)
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


function get_spectrals(; source::Union{Stft, Cwt}, kwargs...)
	_get_spectrals(; s=source.data.spec, sr=source.sr, freq=source.data.freq, kwargs...)
end
