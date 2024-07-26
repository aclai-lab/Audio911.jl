# ---------------------------------------------------------------------------- #
#                             spectral features                                #
# ---------------------------------------------------------------------------- #
mutable struct Spectral
	# setup
    sr::Int64
	freq_range::Tuple{Int64, Int64}
	# data
	centroid::Union{Nothing, AbstractVector{Float64}}
	crest::Union{Nothing, AbstractVector{Float64}}
	decrease::Union{Nothing, AbstractVector{Float64}}
	entropy::Union{Nothing, AbstractVector{Float64}}
	flatness::Union{Nothing, AbstractVector{Float64}}
	flux::Union{Nothing, AbstractVector{Float64}}
	kurtosis::Union{Nothing, AbstractVector{Float64}}
	rolloff::Union{Nothing, AbstractVector{Float64}}
	skewness::Union{Nothing, AbstractVector{Float64}}
	slope::Union{Nothing, AbstractVector{Float64}}
	spread::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function Spectral(;
	sr,
	freq_range = (0, round(Int, sr / 2)),
	centroid = nothing,
	crest = nothing,
	decrease = nothing,
	entropy = nothing,
	flatness = nothing,
	flux = nothing,
	kurtosis = nothing,
	rolloff = nothing,
	skewness = nothing,
	slope = nothing,
	spread = nothing
)
	spectral = Spectral(
        sr,
		freq_range,
		centroid,
		crest,
		decrease,
		entropy,
		flatness,
		flux,
		kurtosis,
		rolloff,
		skewness,
		slope,
		spread
	)
	return spectral
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

function spectral_kurtosis(
	spread::AbstractVector{Float64},
	higher_moment_tmp::AbstractArray{Float64},
	higher_moment_denom::Vector{Float64},
	higher_momement_num::AbstractArray{Float64}
)
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

function spectral_slope(
	s::AbstractArray{Float64},
	arithmetic_mean::AbstractVector{Float64},
	fft_frequencies::AbstractVector{Float64},
)
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
	freq::AbstractVector{Float64},
	spect::Spectral,
)
	# trim to desired frequency range
	s_range = findall(spect.freq_range[1] .<= freq .<= spect.freq_range[2])
	s, freq = s[s_range, :], freq[s_range]

	# common datas
	size_x1 = size(s, 1)
	sum_x1 = vec(sum(s, dims=1))
	arithmetic_mean = sum_x1 ./ size_x1
	spect.centroid = replace!(vec(sum(s .* freq, dims = 1) ./ sum_x1'), NaN => 0)
	higher_moment_tmp = freq .- spect.centroid'
	spect.spread = vec(sqrt.(sum((higher_moment_tmp .^ 2) .* s, dims = 1) ./ sum_x1'))
	higher_moment_denom = (spect.spread .^ 3) .* sum_x1
	higher_momement_num = (higher_moment_tmp .^ 3) .* s

	spect.crest = spectral_crest(s, arithmetic_mean)
	spect.entropy = spectral_entropy(s, sum_x1)
	spect.flatness = spectral_flatness(s, arithmetic_mean)
	spect.flux = spectral_flux(s)
	spect.kurtosis = spectral_kurtosis(spect.spread, higher_moment_tmp, higher_moment_denom, higher_momement_num)
	spect.rolloff = spectral_rolloff(s, freq)
	spect.skewness = spectral_skewness(higher_moment_denom, higher_momement_num)
	spect.decrease = spectral_decrease(s)
	spect.slope = spectral_slope(s, arithmetic_mean, freq)

	return spect
end

function Base.show(io::IO, s::Spectral)
    println(io, "Spectral Features:")
    println(io, "  Sample Rate: $(s.sr) Hz")
    println(io, "  Frequency Range: $(s.freq_range) Hz")
    
    features = [
        ("Centroid", s.centroid),
        ("Crest", s.crest),
        ("Decrease", s.decrease),
        ("Entropy", s.entropy),
        ("Flatness", s.flatness),
        ("Flux", s.flux),
        ("Kurtosis", s.kurtosis),
        ("Rolloff", s.rolloff),
        ("Skewness", s.skewness),
        ("Slope", s.slope),
        ("Spread", s.spread)
    ]
    
    for (name, feature) in features
        status = isnothing(feature) ? "Not computed" : "Computed"
        println(io, "  $name: $status")
    end
end

using Plots

function Base.display(s::Spectral)
    features = [
        ("Centroid", s.centroid),
        ("Crest", s.crest),
        ("Decrease", s.decrease),
        ("Entropy", s.entropy),
        ("Flatness", s.flatness),
        ("Flux", s.flux),
        ("Kurtosis", s.kurtosis),
        ("Rolloff", s.rolloff),
        ("Skewness", s.skewness),
        ("Slope", s.slope),
        ("Spread", s.spread)
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


function get_spectrals(;
    source::Union{Stft, Cwt},
    kwargs...
)
	_get_spectrals(s=source.spec, freq=source.freq, spect=Spectral(; sr=source.sr, kwargs...))
end
