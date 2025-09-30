# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
abstract type AbstractMelSpec end

# ---------------------------------------------------------------------------- #
#                            mel spectrogram struct                            #
# ---------------------------------------------------------------------------- #
struct MelSpec{F,T} <: AbstractMelSpec
	mel_spec :: Matrix{T}
	mel_freq :: Vector{Float64}
	fb       :: Matrix{Float64}
	info     :: NamedTuple

	function MelSpec{F}(
		mel_spec :: Matrix{T},
		mel_freq :: Vector{Float64},
		fb       :: Matrix{Float64},
		info     :: NamedTuple
	) where {F,T}
		new{F,T}(mel_spec, mel_freq, fb, info)
	end
end

# ---------------------------------------------------------------------------- #
#                         mel scale utility functions                          #
# ---------------------------------------------------------------------------- #
function hz2mel(
	hz        :: FreqRange,
	mel_style :: Symbol
)::FreqRange{Float64}
	if mel_style == :htk
		mel = 2595 * log10.(1 .+ reduce(vcat, get_freqs(hz)) / 700)
	else # slaney
		hz = reduce(vcat, get_freqs(hz))
		lin_step = 200 / 3
		log_tep = log(6.4) / 27
		c_point = 1000
		c_point_mel = c_point / lin_step
		is_linear = hz .< c_point
		mel = Float64.(hz)
		mel[is_linear] .= hz[is_linear] / lin_step
		mel[.!is_linear] .= c_point_mel .+ log.(hz[.!is_linear] / c_point) / log_tep
	end
	return FreqRange(mel...)
end

function mel2hz(
	mel       :: LinRange{Float64, Int64},
	mel_style :: Symbol
)
	if mel_style == :htk
		hz = 700 * (exp10.(mel / 2595) .- 1)
	else
		lin_step = 200 / 3
		log_step = log(6.4) / 27
		c_point = 1000
		c_point_mel = c_point / lin_step
		is_linear = mel .< c_point_mel
		hz = [mel;]
		hz[is_linear] .= hz[is_linear] * lin_step
		hz[.!is_linear] .= c_point * exp.(log_step * (mel[.!is_linear] .- c_point_mel))
	end
	return hz
end

# function get_mel_norm_factor(spectrum_type::Symbol, fft_window::Vector{Float64})
# 	if spectrum_type == :power
# 		return 1 / (sum(fft_window)^2)
# 	elseif spectrum_type == :magnitude
# 		return 1 / sum(fft_window)
# 	else
# 		error("Unknown spectrum_type $spectrum_type.")
# 	end
# end

# # ---------------------------------------------------------------------------- #
# #                       semitones scale utility functions                      #
# # ---------------------------------------------------------------------------- #
# # hz2semitone(freq, base_freq) = 12 * log2(freq / base_freq)
# # semitone2hz(z, base_freq)    = base_freq * (2^(z / 12))

# function catch_loudest_index(
# 	lin_spec::AbstractArray{Float64},
# 	freq_array::AbstractVector{Float64},
# 	st_peak_range::Tuple{Int64, Int64},
# )
# 	comp_spec = mean(lin_spec, dims = 1)'

# 	min_index = findfirst((x) -> x >= st_peak_range[1], freq_array)
# 	max_index = findfirst((x) -> x >= st_peak_range[2], freq_array) - 1

# 	loudest_peak, loudest_index = findmax(comp_spec[min_index:max_index])

# 	return loudest_peak, min_index + loudest_index
# end

# ---------------------------------------------------------------------------- #
#                           design fb matrix                           #
# ---------------------------------------------------------------------------- #
### da generalizzare per 
### frequency_scale :mel, :bark, :erb
### filterbanl_design_domain :linear, :warped (da verificare se serve)
function design_fb(
    stft             :: Stft;
    mel_bands        :: Int64=26,
    mel_style        :: Symbol=:htk, 				 # :htk, :slaney, :tuned
    fb_design_domain :: Symbol=:linear,
	fb_norm          :: Symbol=:bandwidth, 		  # :bandwidth, :area, :none
	frequency_scale  :: Symbol=:mel, 			    # TODO :mel, :bark, :erb
    st_peak_range    :: FreqRange=FreqRange(200, 700),
)
    frequency_range = get_info(stft).frequency_range
    stft_size       = get_info(stft).stft_size
    sr              = get_info(stft).sr

	# set the design domain ### da implementare in futuro
	design_domain = fb_design_domain == :linear ? :linear : frequency_scale

	# compute band edges
	# TODO da inserire il caso :erb e :bark

    # quel lin spectrogram andrà studiato bene
	# if mel_style == :tuned
	# 	if isempty(data.lin_spectrogram)
	# 		lin_spectrogram!(setup, data)
	# 	end
	# 	_, loudest_index = catch_loudest_index(data.lin_spectrogram, data.lin_frequencies, st_peak_range)
	# 	melRange = hz2mel((round(Int, data.lin_frequencies[loudest_index]), frequency_range[2]), :htk)
	# else
		melRange = hz2mel(frequency_range, mel_style)
	# end

	# mimic audioflux linear mel_style
	if mel_style == :linear
		lin_fq = collect(0:(stft_size-1)) / stft_size * sr
		band_edges = lin_fq[1:(mel_bands+2)]
	elseif mel_style == :htk || mel_style == :slaney
		band_edges = mel2hz(LinRange(melRange.low, melRange.hi, mel_bands + 2), mel_style)
	# elseif setup.mel_style == :tuned
	# 	band_edges = mel2hz(LinRange(melRange[1], melRange[end], mel_bands + 2), :htk)
	else
		error("Unknown mel_style $(mel_style).")
	end

	### parte esclusiva per mel fb si passa a file designmelfb.m
	# determine the number of bands
	num_edges = length(band_edges)

	# determine the number of valid bands
	valid_num_edges = sum((band_edges .- (sr / 2)) .< sqrt(eps(Float64)))
	valid_num_bands = valid_num_edges - 2

	# preallocate the filter bank
	fb = zeros(Float64, stft_size, mel_bands)
	mel_freq = band_edges[2:(end-1)]

	# Set this flag to true if the number of FFT length is insufficient to
	# compute the specified number of mel bands
	FFTLengthTooSmall = false

	# if :hz 
	linFq = collect(0:(stft_size-1)) / stft_size * sr

	# Determine inflection points
	@assert(valid_num_edges <= num_edges)
	p = zeros(Float64, valid_num_edges, 1)

	for edge_n in 1:valid_num_edges
		for index in eachindex(linFq)
			if linFq[index] > band_edges[edge_n]
				p[edge_n] = index
				break
			end
		end
	end

	FqMod = linFq

	# Create triangular filters for each band
	bw = diff(band_edges)

	for k in 1:Int(valid_num_bands)
		# Rising side of triangle
		for j in Int(p[k]):(Int(p[k+1])-1)
			fb[j, k] = (FqMod[j] - band_edges[k]) / bw[k]
		end
		# Falling side of triangle
		for j in Int(p[k+1]):(Int(p[k+2])-1)
			fb[j, k] = (band_edges[k+2] - FqMod[j]) / bw[k+1]
		end
		emptyRange1 = p[k] .> p[k+1] - 1
		emptyRange2 = p[k+1] .> p[k+2] - 1
		if (!FFTLengthTooSmall && (emptyRange1 || emptyRange2))
			FFTLengthTooSmall = true
		end
	end

	# mirror two sided
	range = get_onesided_fft_range(stft_size)
	range = range[2:end]
	fb[end:-1:(end-length(range)+1), :] = fb[range, :]

	fb = fb'

	# normalizzazione    
	BW = band_edges[3:end] - band_edges[1:(end-2)]

	if (fb_norm == :area)
		weight_per_band = sum(fb, dims = 2)
		if frequency_scale != :erb
			weight_per_band = weight_per_band / 2
		end
	elseif (fb_norm == :bandwidth)
		weight_per_band = BW / 2
	else
		weight_per_band = ones(1, mel_bands)
	end

	for i in 1:(mel_bands)
		if (weight_per_band[i] != 0)
			fb[i, :] = fb[i, :] ./ weight_per_band[i]
		end
	end

	# get one side
	range = get_onesided_fft_range(stft_size)
	fb = fb[:, range]
	# manca la parte relativa a :erb e :bark

	# setta fattore di normalizzazione
	# if setup.window_norm
	# 	win_norm_factor = get_mel_norm_factor(setup.spectrum_type, data.fft_window)
	# 	fb = fb * win_norm_factor
	# end

	return fb, mel_freq
end

# ---------------------------------------------------------------------------- #
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
function get_melspec(
	stft::Stft;
    mel_bands        :: Int64=26,
	mel_style        :: Symbol=:htk, 				 # :htk, :slaney, :tuned
	fb_design_domain :: Symbol=:linear,
	fb_norm          :: Symbol=:bandwidth, 		  # :bandwidth, :area, :none
	frequency_scale  :: Symbol=:mel, 			    # TODO :mel, :bark, :erb
	st_peak_range    :: FreqRange=FreqRange(200, 700),
)
    data = get_stft(stft)
    win_size = get_info(stft).win_size
    win_step = get_info(stft).win_step

	fb, mel_freq = design_fb(
        stft;
        mel_bands,
        mel_style,
        fb_design_domain,
        fb_norm,
        frequency_scale,
        st_peak_range
    )

	num_hops = size(data, 2)

	# apply fb
	# if (setup.spectrum_type == :power)
	mel_spec = reshape(fb * data, mel_bands, num_hops)
	# else
	#     #TODO
	#     error("magnitude not yet implemented.")
	# end

	# mel_spec = transpose(mel_spec)

    return mel_spec, fb, mel_freq
end

function get_melspec(
	frames          :: AudioFrames;
	stft_size       :: Int64=get_wsize(frames),
	frequency_range :: FreqRange=FreqRange(0, get_info(frames).sr÷2),
	spectrum_type   :: Symbol=:power, # :power, :magnitude
    kwargs...
)
    stft = get_stft(frames; stft_size, frequency_range, spectrum_type)
    get_melspec(stft; kwargs...)
end

function get_melspec(
	afile :: AudioFile;
    win   :: WinFunction=MovingWindow(
                            window_size=samplerate(afile)≤8000 ? 256 : 512,
                            window_step=samplerate(afile)≤8000 ? 128 : 256
                        ),
	type  :: Tuple{Symbol, Symbol}=(:hann, :periodic),
	kwargs...
	)
	frames = get_frames(afile; win, type)
	get_melspec(frames; kwargs...)
end

# # ---------------------------------------------------------------------------- #
# #                          logaritmic mel spectrogram                          #
# # ---------------------------------------------------------------------------- #
# # TODO prova a fare le delta del log mel
# function get_log_mel!(
# 	setup::AudioSetup,
# 	data::AudioData,
# )
# 	# Reference:
# 	# https://dsp.stackexchange.com/questions/85501/log-of-fb-energies

# 	# Rectify
# 	mel_spec = deepcopy(data.mel_spectrogram')

# 	if setup.normalization_type == :standard
# 		mel_spec[mel_spec.==0] .= floatmin(Float64)
# 	elseif setup.normalization_type == :dithered
# 		mel_spec[mel_spec.<1e-8] .= 1e-8
# 	else
# 		@warn("Unknown $setup.normalization_type normalization type, defaulting to standard.")
# 		mel_spec[mel_spec.==0] .= floatmin(Float64)
# 	end

# 	data.log_mel = log10.(mel_spec)'
# end

# # ---------------------------------------------------------------------------- #
# #                             mfcc related functions                           #
# # ---------------------------------------------------------------------------- #
# function create_DCT_matrix(
# 	mel_coeffs::Int64,
# )
# 	# create DCT matrix
# 	matrix = zeros(Float64, mel_coeffs, mel_coeffs)
# 	s0 = sqrt(1 / mel_coeffs)
# 	s1 = sqrt(2 / mel_coeffs)
# 	piCCast = 2 * pi / (2 * mel_coeffs)

# 	matrix[1, :] .= s0
# 	for k in 1:mel_coeffs, n in 2:mel_coeffs
# 		matrix[n, k] = s1 * cos(piCCast * (n - 1) * (k - 0.5))
# 	end

# 	matrix
# end

# function audioDelta(
# 	x::AbstractMatrix{T},
# 	window_length::Int64,
# 	source::Symbol = :standard,
# ) where {T <: AbstractFloat}

# 	# define window shape
# 	m = Int(floor(window_length / 2))
# 	b = collect(m:-1:(-m)) ./ sum((1:m) .^ 2)

# 	if source == :transposed
# 		filt(b, 1.0, x')'   #:audioflux setting
# 	else
# 		filt(b, 1.0, x)     #:matlab setting
# 	end
# end

# function get_mfcc!(
# 	setup::AudioSetup,
# 	data::AudioData,
# )
# 	# Rectify
# 	mel_spec = deepcopy(data.mel_spectrogram')

# 	if setup.normalization_type == :standard
# 		mel_spec[mel_spec.==0] .= floatmin(Float64)
# 	elseif setup.normalization_type == :dithered
# 		mel_spec[mel_spec.<1e-8] .= 1e-8
# 	else
# 		@warn("Unknown $setup.normalization_type normalization type, defaulting to standard.")
# 		mel_spec[mel_spec.==0] .= floatmin(Float64)
# 	end

# 	# Design DCT matrix
# 	DCTmatrix = create_DCT_matrix(setup.mel_bands)

# 	# apply DCT matrix
# 	if (setup.rectification == :log)
# 		coeffs = DCTmatrix * log10.(mel_spec)
# 	elseif (setup.rectification == :cubic_root)
# 		# apply DCT matrix
# 		coeffs = DCTmatrix * mel_spec .^ (1 / 3)
# 	else
# 		@warn("Unknown $rectification DCT matrix rectification, defaulting to log.")
# 		coeffs = DCTmatrix * log10.(mel_spec)
# 	end

# 	# reduce to mfcc coefficients
# 	data.mfcc_coeffs = coeffs[1:(setup.num_coeffs), :]'

# 	# log energy calc
# 	if setup.log_energy_source == :mfcc
# 		log_energy = sum(eachrow(mel_spec .^ 2)) / setup.mel_bands

# 		if setup.normalization_type == :standard
# 			log_energy[log_energy.==0] .= floatmin(Float64)
# 		elseif setup.normalization_type == :dithered
# 			log_energy[log_energy.<1e-8] .= 1e-8
# 		end

# 		data.log_energy = log.(log_energy)
# 	end

# 	if (setup.log_energy_pos == :append)
# 		data.mfcc_coeffs = hcat(data.mfcc_coeffs, data.log_energy)
# 	elseif (setup.log_energy_pos == :replace)
# 		data.mfcc_coeffs = hcat(data.log_energy, data.mfcc_coeffs[:, 2:end])
# 	end
# end

# function get_mfcc_deltas!(
# 	setup::AudioSetup,
# 	data::AudioData,
# )
# 	data.mfcc_delta = audioDelta(
# 		data.mfcc_coeffs, setup.delta_window_length, setup.delta_matrix)
# 	data.mfcc_deltadelta = audioDelta(
# 		data.mfcc_delta, setup.delta_window_length, setup.delta_matrix)
# end
