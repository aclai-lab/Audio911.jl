# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct MelSpecSetup <: AbstractSetup
    sr            :: Int64
    freq_range    :: FreqRange
    spectrum_type :: Base.Callable
	win_norm      :: Bool
end

# ---------------------------------------------------------------------------- #
#                          linear spectrogram struct                           #
# ---------------------------------------------------------------------------- #
struct MelSpec{F,T} <: AbstractSpectrogram
    spec :: Matrix{T}
    freq :: StepRangeLen
    info :: MelSpecSetup

    function MelSpec{F}(
        spec :: Matrix{T},
        freq :: StepRangeLen,
        info :: MelSpecSetup
    ) where {F,T}
        new{F,T}(spec, freq, info)
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
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
function LinSpec(
	stft       :: Stft;
	freq_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr>>1),
	win_norm   :: Bool=false
)::LinSpec
	spec = get_spec(stft)
	freq = get_freq(stft)

	sr            = get_sr(stft)
	stft_size     = get_stft_size(stft)
	spectrum_type = get_spec_type(stft)
	window        = get_window(stft)

	freq_range isa Tuple && (freq_range = FreqRange(first(freq_range), last(freq_range)))

	if freq_range != FreqRange(0, sr >> 1)
		bin_low, bin_high = get_freq_range(freq_range, stft_size, sr)
		spec = @views spec[bin_low:bin_high, :]
		freq = freq[bin_low:bin_high]
	end

	win_norm_func = eval(Symbol("win" * string(spectrum_type)))
	win_norm && (spec = win_norm_func(spec, window))

	info = LinSpecInfo(sr, freq_range, spectrum_type, win_norm)

	return LinSpec{typeof(stft)}(spec .* 2, freq, info)
end

function get_melspec(
	stft::Stft;
	freq_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr>>1),
	win_norm   :: Bool=false,
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
