#------------------------------------------------------------------------------#
#                              define structures                               #
#------------------------------------------------------------------------------#
function audio_setup(
	sr::Int64;

	# stft
	stft_length::Int64 = sr <= 8000 ? 256 : 512,
	win::AbstractVector{Float64} = Float64[],
	win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
	win_length::Int64 = stft_length, 					# standard setting: round(Int, 0.03 * sr)
	overlap_length::Int64 = round(Int, stft_length / 2), # standard setting: round(Int, 0.02 * sr)
	spec_norm::Symbol = :power, # :none, :power, :magnitude, :pow2mag

	# lin
	win_norm::Symbol = :none, # :none, :power, :magnitude
	freq_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
	apply_log::Bool = false,

	# mel
	mel_style::Symbol = :htk, 							# :htk, :slaney, :tuned
	mel_bands::Int64 = 26,
	filterbank_design_domain::Symbol = :linear,
	filterbank_normalization::Symbol = :bandwidth, 		# :bandwidth, :area, :none
	frequency_scale::Symbol = :mel, 					# TODO :mel, :bark, :erb
	st_peak_range::Tuple{Int64, Int64} = (200, 700),

	# chroma
	bins_octave::Int64 = 12,
	center_freq::Int64 = 1000,
	gaussian_sd::Int64 = 1,

	# mfcc
	mfcc_coeffs::Int64 = 13,
	normalization_type::Symbol = :dithered, 			# :standard, :dithered
	rectification::Symbol = :log, 						# :log, :cubic_root
	log_energy_source::Symbol = :standard, 				# :standard (after windowing), :mfcc
	log_energy_pos::Symbol = :none, 					#:append, :replace, :none
	delta_win_length::Int64 = 9,
	delta_matrix::Symbol = :transposed, 				# :standard, :transposed

	# spectral
	spectral_spectrum::Symbol = :lin, 					# :lin, :mel

	# f0
	f0_method::Symbol = :nfc,
	f0_range::Tuple{Int64, Int64} = (50, 400),
	median_filter_length::Int64 = 1,

	# constant-q transform	
	freq_limits::Tuple{Float64, Float64} = (0.0, 0.0),
	transform_type::Symbol = :full,
)
	# TODO metti warning ed errori

	AudioSetup(
		sr = sr,

		stft = StftSetup(
			stft_length = stft_length,
			win = win,
			win_type = win_type,
			win_length = win_length,
			overlap_length = overlap_length,
			spec_norm = spec_norm,
		),

		lin_spec = LinSetup(
			win_norm = win_norm,
			freq_range = freq_range,
			apply_log = apply_log
		),

		# mel
		mel_style = mel_style,
		mel_bands = mel_bands,
		filterbank_design_domain = filterbank_design_domain,
		filterbank_normalization = filterbank_normalization,
		frequency_scale = frequency_scale,
		st_peak_range = st_peak_range,

		# chroma
		bins_octave = bins_octave,
		center_freq = center_freq,
		gaussian_sd = gaussian_sd,

		# mfcc
		mfcc_coeffs = mfcc_coeffs,
		normalization_type = normalization_type,
		rectification = rectification,
		log_energy_source = log_energy_source,
		log_energy_pos = log_energy_pos,
		delta_win_length = delta_win_length,
		delta_matrix = delta_matrix,

		# spectral
		spectral_spectrum = spectral_spectrum,

		# f0
		f0_method = f0_method,
		f0_range = f0_range,
		median_filter_length = median_filter_length,

		# constant-q transform
		freq_limits = freq_limits,
		transform_type = transform_type,
	)
end

#------------------------------------------------------------------------------#
#                                audio object                                  #
#------------------------------------------------------------------------------#

function audio_obj(
	x::AbstractVector{Float64},
	sr::Int64;
	preemphasis::Union{Float64, Nothing} = nothing,
	kwargs...,
)
	setup = audio_setup(sr; kwargs...)

	if preemphasis !== nothing
		zi = 2 * x[1] - x[2]
		filt!(x, [1.0, - preemphasis], 1.0, x, [zi])
		# # aclai preemphasis
		# x = filt(PolynomialRatio([1.0, -preemphasis], [1.0]), x)
	end

	# check x audio sample length
	if length(x) == 0
		@warn("Sample is empty, skipped.")
		return nothing
	elseif length(x) < setup.stft.stft_length
		@warn("length of sample is < than stft windows, skipped.")
		return nothing
	else
		data = AudioData(x = x)
		audio_obj = AudioObj(setup, data)
	end

	return audio_obj
end

function audio_obj(
	x::AbstractVector{T},
	sr::Int64;
	preemphasis = nothing,
	kwargs...,
) where {T <: AbstractFloat}
	audio_obj(Float64.(x), sr; preemphasis, kwargs...)
end

function audio_obj(
	filepath::String,
	sr::Int64;
	preemphasis = nothing,
	kwargs...,
)
	x, sr = load_audio(filepath, sr)
	x = normalize_audio(x)
	audio_obj(Float64.(x), sr; preemphasis, kwargs...)
end

#------------------------------------------------------------------------------#
#                            get features functions                            #
#------------------------------------------------------------------------------#
function get_stft(audio_obj::AudioObj)
	get_stft!(audio_obj)
	return audio_obj.data.stft.stft
end

function get_mel_spec(audio_obj::AudioObj)
	if isempty(audio_obj.data.mel_spectrogram)
		if isempty(audio_obj.data.stft.stft)
			get_stft!(audio_obj)
		end
		get_mel_spec!(audio_obj.setup, audio_obj.data)
	end
	return audio_obj.data.mel_spectrogram, audio_obj.data.mel_frequencies
end

function get_log_mel(audio_obj::AudioObj)
	if isempty(audio_obj.data.log_mel)
		if isempty(audio_obj.data.mel_spectrogram)
			if isempty(audio_obj.data.stft.stft)
				get_stft!(audio_obj)
			end
			get_mel_spec(audio_obj)
		end
		get_log_mel!(audio_obj.setup, audio_obj.data)
	end
	return audio_obj.data.log_mel
end

function get_mfcc(audio_obj::AudioObj)
	if isempty(audio_obj.data.mfcc_coeffs)
		if isempty(audio_obj.data.mel_spectrogram)
			if isempty(audio_obj.data.stft.stft)
				get_stft!(audio_obj)
			end
			get_mel_spec(audio_obj)
		end
		get_mfcc!(audio_obj.setup, audio_obj.data)
	end
	return audio_obj.data.mfcc_coeffs
end

function get_mfcc_delta(audio_obj::AudioObj)
	if isempty(audio_obj.data.mfcc_delta)
		if isempty(audio_obj.data.mfcc_coeffs)
			if isempty(audio_obj.data.mel_spectrogram)
				if isempty(audio_obj.data.stft.stft)
					get_stft!(audio_obj)
				end
				get_mel_spec(audio_obj)
			end
			get_mfcc(audio_obj)

		end
		get_mfcc_deltas!(audio_obj.setup, audio_obj.data)
	end

	return hcat((audio_obj.data.mfcc_delta, audio_obj.data.mfcc_deltadelta)...)
end

function get_spectrals(audio_obj::AudioObj)
	if audio_obj.setup.spectral_spectrum == :lin && isempty(audio_obj.data.stft.stft)
		get_stft!(audio_obj)

	elseif audio_obj.setup.spectral_spectrum == :mel && isempty(audio_obj.data.mel_spectrogram)
		if isempty(audio_obj.data.stft.stft)
			get_stft!(audio_obj)
		end
		mel_spectrogram!(audio_obj)
	end
	if isempty(audio_obj.data.spectral_centroid)
		get_spectrals!(audio_obj.setup, audio_obj.data)
	end

	return hcat(
		(
			audio_obj.data.spectral_centroid,
			audio_obj.data.spectral_crest,
			audio_obj.data.spectral_decrease,
			audio_obj.data.spectral_entropy,
			audio_obj.data.spectral_flatness,
			audio_obj.data.spectral_flux,
			audio_obj.data.spectral_kurtosis,
			audio_obj.data.spectral_rolloff,
			audio_obj.data.spectral_skewness,
			audio_obj.data.spectral_slope,
			audio_obj.data.spectral_spread,
		)...,
	)
end

function get_f0(audio_obj::AudioObj)
	if isempty(audio_obj.data.f0)
		get_f0!(audio_obj.setup, audio_obj.data)
	end

	return audio_obj.data.f0
end

function get_cqt(audio_obj::AudioObj)
	if isempty(audio_obj.data.cqt_spec)
		get_cqt_spec!(audio_obj.setup, audio_obj.data)
	end

	return audio_obj.data.cqt_spec
end

function get_full(audio_obj::AudioObj)
	get_stft!(audio_obj)
	get_mel_spec(audio_obj)
	# get_log_mel(audio_obj)
	get_mfcc(audio_obj)
	get_mfcc_delta(audio_obj)
	get_spectrals(audio_obj)
	get_f0(audio_obj)

	return hcat(
		(
			audio_obj.data.mel_spectrogram,
			# audio_obj.data.log_mel,
			audio_obj.data.mfcc_coeffs,
			audio_obj.data.mfcc_delta,
			audio_obj.data.mfcc_deltadelta,
			audio_obj.data.spectral_centroid,
			audio_obj.data.spectral_crest,
			audio_obj.data.spectral_decrease,
			audio_obj.data.spectral_entropy,
			audio_obj.data.spectral_flatness,
			audio_obj.data.spectral_flux,
			audio_obj.data.spectral_kurtosis,
			audio_obj.data.spectral_rolloff,
			audio_obj.data.spectral_skewness,
			audio_obj.data.spectral_slope,
			audio_obj.data.spectral_spread,
			audio_obj.data.f0,
		)...)
end

#------------------------------------------------------------------------------#
#                               Audio911 caller                                #
#------------------------------------------------------------------------------#
function get_features(
	audio_obj::Union{AudioObj, Nothing},
	feat::Symbol = :full,
)
	func_call = Dict([
		:full => get_full, 
		:stft => get_stft,
		:mel => get_mel_spec,
		:logmel => get_log_mel,
		:mfcc => get_mfcc,
		:mfcc_delta => get_mfcc_delta,
		:spectral => get_spectrals,
		:f0 => get_f0,
		:cqt => get_cqt,

		:gender_set =>
		(x) -> begin
			get_mfcc(x)
			get_f0(x)
			return hcat(
				(
					audio_obj.data.mfcc_coeffs,
					audio_obj.data.f0,
				)...)
		end,
		
		:age_set =>
			(x) -> begin
				get_full(x)
				return hcat(
					(
						# audio_obj.data.mel_spectrogram,
						# audio_obj.data.log_mel,
						audio_obj.data.mfcc_coeffs,
						audio_obj.data.mfcc_delta,
						audio_obj.data.mfcc_deltadelta,
						# audio_obj.data.spectral_centroid,
						# audio_obj.data.spectral_crest,
						# audio_obj.data.spectral_decrease,
						# audio_obj.data.spectral_entropy,
						# audio_obj.data.spectral_flatness,
						# audio_obj.data.spectral_flux,
						# audio_obj.data.spectral_kurtosis,
						# audio_obj.data.spectral_rolloff,
						# audio_obj.data.spectral_skewness,
						# audio_obj.data.spectral_slope,
						# audio_obj.data.spectral_spread,
						# audio_obj.data.f0,
					)...)
			end])

	if !isnothing(audio_obj)
		if haskey(func_call, feat)
			func_call[feat](audio_obj)
		else
			println("feature $feat not available.")
			return nothing
		end
	else
		return nothing
	end
end

# function get_features(
# 	x::AbstractVector{T},
# 	sr::Int64,
# 	feat::Symbol = :full;
# 	kwargs...,
# ) where {T <: AbstractFloat}
# 	get_features(audio_obj(x, sr; kwargs...), feat)
# end

# function get_features(
# 	filepath::String,
# 	sr::Int64,
# 	feat::Symbol = :full;
# 	kwargs...,
# )
# 	get_features(audio_obj(filepath, sr; kwargs...), feat)
# end