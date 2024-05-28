"""
Audio911

può essere utilizzata in 2 modi differenti:
#########################################################################################################
1-creare un oggetto audio_obj

audio = audio_features_obj(x, sr)

audio.get_fft()
audio.get_lin_spec()
audio.get_mel_spec()
audio.get_mfcc()
audio.get_spectrals()
audio.get_f0()
audio.get_features(profile, per ora solo :full)

#########################################################################################################
2-utilizzare l'invocazione get_feature per ottenere le features separatamente

costruzione:
get_features(
	x AbstractFloat file audio mono
	sr Int64 frequenza di campionamento
	feat Symbol audio feature da estrarre
		:fft fast fourier transform (da implementare meglio)
		:lin spettrogramma lineare
		:mel spettrogramma mel
		:mfcc coefficienti mfcc e relative delta e deltadelta
		:spectrals le feature spettrali:
		centroid, crest, decrease, entropy, flatness, flux, kurtosis, rolloff, skewness, slope, spread
		:f0 frequenza fondamentale
	kwargs...

##############################################################################################

da spostare in signal data structure 

paramentri addizionali
sia per l'oggetto audio, che per la chiamata a feature singola

# fft
stft_length::Int64 = 256,
dimensione finestra fft, valori consigliati: 256, 512, 1024

win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
win_length::Int64 = stft_length,
overlap_length::Int64 = round(Int, stft_length * 0.500),
parametri relativi alla finestrazione della fft
di default audio911 usa una finestra di tipo hann, anzichè hamming
con una finestra della stessa dimensione della finestra fft
e un overlap pari alla metà del suo valore

i valori standard sarebbero questi
# win_length::Int64 = Int(round(0.03 * sr)),
# overlap_length::Int64 = Int(round(0.02 * sr)),

window_norm::Bool = false, normalizzazione delle finestre

# spectrum
freq_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
limiti banda, importantissimi per isolare la porzione di spettro dove si prevede di recuperare l'informazione

spectrum_type::Symbol = :power, # :power, :magnitude
tipo di spettro, di default :power, molto raramente si usa magnitude

# mel
mel_style::Symbol = :htk, # :htk, :slaney
tipo di banco filtro dello spettrogramma mel
la tipologia htk surclassa la tipologia slaney in tutti i nostri esperimenti

mel_bands::Int64 = 26,
numero di bande che compongono lo spettrogramma mel. 26 è il valore di default
ma non c'è un vero e proprio standard di questo valore.
da ricordare che nel caso si voglia utilizzare anche la mfcc,
i coefficienti della mfcc vengono calcolati sulle prime bande dello spettrogramma
quindi una variazione di questo valore comporta anche un diverso funzionamento della mfcc
il cui valore "num_coeff" andrà opportunamente tarato.

filterbank_design_domain::Symbol = :linear,
filterbank_normalization::Symbol = :bandwidth, # :bandwidth, :area, :none
frequency_scale::Symbol = :mel,
l'implementazione corretta di questi paramentri è da completare

# mfcc
mfcc_coeffs::Int64 = 13,
numero delle bande mfcc, vedi sopra

normalization_type::Symbol = :dithered, # :standard, :dithered
paramentro preso da audioflux:
i valori dei coefficienti mfcc, se vanno sotto la soglia di 1e-8
vengono normalizzati a questo valore.
mentre matlab non ha una soglia limite.
dithered > audiofluz, standard > matlab

rectification::Symbol = :log,
log_energy_source::Symbol = :standard, # :standard (after windowing), :mfcc
log_energy_pos::Symbol = :none, #:append, :replace, :none
spesso viene salvato, nella mfcc, il valore del volume in log.
questi parametri definiscono dove viene calcolata e dove salvarla: 
:append viene creato un n-esimo coefficiente, :replace la log energy va a sostituire il primo coefficiente mfcc
se ne sconsiglia comunque l'uso.

delta_win_length::Int64 = 9,
finestra di calcolo della derivata

delta_matrix::Symbol = :transposed, # :standard, :transposed
preso da audioflux che calcola le delta sull'asse delle frequenze anzichè sull'asse temporale
potrebbe sembrare un errore, ma potrebbe anche non esserlo

# spectral
spectral_spectrum::Symbol = :lin, # :lin, :mel
si può scegliere su che spettrogramma calcolare le spectral features: se partendo dal lineare o dal mel

# f0
f0_method::Symbol = :nfc,
f0_range::Tuple{Int64, Int64} = (50, 400),
median_filter_length::Int64 = 1
Questi paramentri sono in fase di studio
"""

#------------------------------------------------------------------------------#
#                              define structures                               #
#------------------------------------------------------------------------------#
function audio_setupdev(
	sr::Int64;

	# stft
	stft_length::Int64 = sr <= 8000 ? 256 : 512,
	win::AbstractVector{Float64} = Float64[],
	win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
	win_length::Int64 = stft_length, 					# standard setting: round(Int, 0.03 * sr)
	overlap_length::Int64 = round(Int, stft_length / 2), # standard setting: round(Int, 0.02 * sr)
	spec_norm::Symbol = :power, #:power, :magnitude, :winpower, :winmagnitude
	freq_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),


	# # mel
	# mel_style::Symbol = :htk, 							# :htk, :slaney, :tuned
	# mel_bands::Int64 = 26,
	# filterbank_design_domain::Symbol = :linear,
	# filterbank_normalization::Symbol = :bandwidth, 		# :bandwidth, :area, :none
	# frequency_scale::Symbol = :mel, 					# TODO :mel, :bark, :erb
	# st_peak_range::Tuple{Int64, Int64} = (200, 700),

	# # chroma
	# bins_octave::Int64 = 12,
	# center_freq::Int64 = 1000,
	# gaussian_sd::Int64 = 1,

	# # mfcc
	# mfcc_coeffs::Int64 = 13,
	# normalization_type::Symbol = :dithered, 			# :standard, :dithered
	# rectification::Symbol = :log, 						# :log, :cubic_root
	# log_energy_source::Symbol = :standard, 				# :standard (after windowing), :mfcc
	# log_energy_pos::Symbol = :none, 					#:append, :replace, :none
	# delta_win_length::Int64 = 9,
	# delta_matrix::Symbol = :transposed, 				# :standard, :transposed

	# # spectral
	# spectral_spectrum::Symbol = :lin, 					# :lin, :mel

	# # f0
	# f0_method::Symbol = :nfc,
	# f0_range::Tuple{Int64, Int64} = (50, 400),
	# median_filter_length::Int64 = 1,

	# # constant-q transform	
	# freq_limits::Tuple{Float64, Float64} = (0.0, 0.0),
	# transform_type::Symbol = :full,
)
	# TODO metti warning ed errori

	AudioSetupDev(
		sr = sr,

		stft = StftSetup(
			stft_length = stft_length,
			win = win,
			win_type = win_type,
			win_length = win_length,
			overlap_length = overlap_length,
			spec_norm = spec_norm,
			freq_range = freq_range,
		)

		# # mel
		# mel_style = mel_style,
		# mel_bands = mel_bands,
		# filterbank_design_domain = filterbank_design_domain,
		# filterbank_normalization = filterbank_normalization,
		# frequency_scale = frequency_scale,
		# st_peak_range = st_peak_range,

		# # chroma
		# bins_octave = bins_octave,
		# center_freq = center_freq,
		# gaussian_sd = gaussian_sd,

		# # mfcc
		# mfcc_coeffs = mfcc_coeffs,
		# normalization_type = normalization_type,
		# rectification = rectification,
		# log_energy_source = log_energy_source,
		# log_energy_pos = log_energy_pos,
		# delta_win_length = delta_win_length,
		# delta_matrix = delta_matrix,

		# # spectral
		# spectral_spectrum = spectral_spectrum,

		# # f0
		# f0_method = f0_method,
		# f0_range = f0_range,
		# median_filter_length = median_filter_length,

		# # constant-q transform
		# freq_limits = freq_limits,
		# transform_type = transform_type,
	)
end

#------------------------------------------------------------------------------#
#                                audio object                                  #
#------------------------------------------------------------------------------#

function audio_objdev(
	x::AbstractVector{Float64},
	sr::Int64;
	preemphasis::Union{Float64, Nothing} = nothing,
	kwargs...,
)
	setup = audio_setupdev(sr; kwargs...)

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
		data = AudioDataDev(x = x)
		audio_obj = AudioObjDev(setup, data)
	end

	return audio_obj
end

function audio_objdev(
	x::AbstractVector{T},
	sr::Int64;
	preemphasis = nothing,
	kwargs...,
) where {T <: AbstractFloat}
	audio_objdev(Float64.(x), sr; preemphasis, kwargs...)
end

function audio_objdev(
	filepath::String,
	sr::Int64;
	preemphasis = nothing,
	kwargs...,
)
	x, sr = load_audio(filepath, sr)
	x = normalize_audio(x)
	audio_objdev(Float64.(x), sr; preemphasis, kwargs...)
end

# #------------------------------------------------------------------------------#
# #                            get features functions                            #
# #------------------------------------------------------------------------------#
# function get_fft(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.fft)
# 		get_fft!(audio_obj.setup, audio_obj.data)
# 	end
# 	return audio_obj.data.fft'
# end

# function get_lin_spec(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.lin_spectrogram)
# 		if isempty(audio_obj.data.fft)
# 			get_fft(audio_obj)
# 		end
# 		lin_spectrogram!(audio_obj.setup, audio_obj.data)
# 	end
# 	return audio_obj.data.lin_spectrogram, audio_obj.data.lin_frequencies
# end

# function get_mel_spec(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.lin_spectrogram)
# 		if isempty(audio_obj.data.fft)
# 			get_fft(audio_obj)
# 		end
# 		get_mel_spec!(audio_obj.setup, audio_obj.data)
# 	end
# 	return audio_obj.data.mel_spectrogram, audio_obj.data.mel_frequencies
# end

# function get_log_mel(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.log_mel)
# 		if isempty(audio_obj.data.mel_spectrogram)
# 			if isempty(audio_obj.data.fft)
# 				get_fft(audio_obj)
# 			end
# 			get_mel_spec(audio_obj)
# 		end
# 		get_log_mel!(audio_obj.setup, audio_obj.data)
# 	end
# 	return audio_obj.data.log_mel
# end

# function get_mfcc(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.mfcc_coeffs)
# 		if isempty(audio_obj.data.mel_spectrogram)
# 			if isempty(audio_obj.data.fft)
# 				get_fft(audio_obj)
# 			end
# 			get_mel_spec(audio_obj)
# 		end
# 		get_mfcc!(audio_obj.setup, audio_obj.data)
# 	end
# 	return audio_obj.data.mfcc_coeffs
# end

# function get_mfcc_delta(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.mfcc_delta)
# 		if isempty(audio_obj.data.mfcc_coeffs)
# 			if isempty(audio_obj.data.mel_spectrogram)
# 				if isempty(audio_obj.data.fft)
# 					get_fft(audio_obj)
# 				end
# 				get_mel_spec(audio_obj)
# 			end
# 			get_mfcc(audio_obj)

# 		end
# 		get_mfcc_deltas!(audio_obj.setup, audio_obj.data)
# 	end

# 	return hcat((audio_obj.data.mfcc_delta, audio_obj.data.mfcc_deltadelta)...)
# end

# function get_spectrals(audio_obj::AudioObjDev)
# 	if audio_obj.setup.spectral_spectrum == :lin && isempty(audio_obj.data.lin_spectrogram)
# 		if isempty(audio_obj.data.fft)
# 			get_fft(audio_obj)
# 		end
# 		get_lin_spec(audio_obj)
# 	elseif audio_obj.setup.spectral_spectrum == :mel && isempty(audio_obj.data.mel_spectrogram)
# 		if isempty(audio_obj.data.fft)
# 			get_fft(audio_obj)
# 		end
# 		mel_spectrogram!(audio_obj)
# 	end
# 	if isempty(audio_obj.data.spectral_centroid)
# 		get_spectrals!(audio_obj.setup, audio_obj.data)
# 	end

# 	return hcat(
# 		(
# 			audio_obj.data.spectral_centroid,
# 			audio_obj.data.spectral_crest,
# 			audio_obj.data.spectral_decrease,
# 			audio_obj.data.spectral_entropy,
# 			audio_obj.data.spectral_flatness,
# 			audio_obj.data.spectral_flux,
# 			audio_obj.data.spectral_kurtosis,
# 			audio_obj.data.spectral_rolloff,
# 			audio_obj.data.spectral_skewness,
# 			audio_obj.data.spectral_slope,
# 			audio_obj.data.spectral_spread,
# 		)...,
# 	)
# end

# function get_f0(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.f0)
# 		get_f0!(audio_obj.setup, audio_obj.data)
# 	end

# 	return audio_obj.data.f0
# end

# function get_cqt(audio_obj::AudioObjDev)
# 	if isempty(audio_obj.data.cqt_spec)
# 		get_cqt_spec!(audio_obj.setup, audio_obj.data)
# 	end

# 	return audio_obj.data.cqt_spec
# end

# function get_full(audio_obj::AudioObjDev)
# 	get_fft(audio_obj)
# 	get_mel_spec(audio_obj)
# 	get_log_mel(audio_obj)
# 	get_mfcc(audio_obj)
# 	get_mfcc_delta(audio_obj)
# 	if audio_obj.setup.spectral_spectrum == :lin
# 		get_lin_spec(audio_obj)
# 	end
# 	get_spectrals(audio_obj)
# 	get_f0(audio_obj)

# 	return hcat(
# 		(
# 			audio_obj.data.mel_spectrogram,
# 			# audio_obj.data.log_mel,
# 			audio_obj.data.mfcc_coeffs,
# 			audio_obj.data.mfcc_delta,
# 			audio_obj.data.mfcc_deltadelta,
# 			audio_obj.data.spectral_centroid,
# 			audio_obj.data.spectral_crest,
# 			audio_obj.data.spectral_decrease,
# 			audio_obj.data.spectral_entropy,
# 			audio_obj.data.spectral_flatness,
# 			audio_obj.data.spectral_flux,
# 			audio_obj.data.spectral_kurtosis,
# 			audio_obj.data.spectral_rolloff,
# 			audio_obj.data.spectral_skewness,
# 			audio_obj.data.spectral_slope,
# 			audio_obj.data.spectral_spread,
# 			audio_obj.data.f0,
# 		)...)
# end

#------------------------------------------------------------------------------#
#                               Audio911 caller                                #
#------------------------------------------------------------------------------#
# function get_features(
# 	audio_obj::Union{AudioObj, Nothing},
# 	feat::Symbol = :full,
# )
# 	func_call = Dict([
# 		:full => get_full, 
# 		:fft => get_fft,
# 		:lin => get_lin_spec,
# 		:mel => get_mel_spec,
# 		:logmel => get_log_mel,
# 		:mfcc => get_mfcc,
# 		:mfcc_delta => get_mfcc_delta,
# 		:spectral => get_spectrals,
# 		:f0 => get_f0,
# 		:cqt => get_cqt,
# 		:gender_set =>
# 		(x) -> begin
# 			get_mfcc(x)
# 			get_f0(x)
# 			return hcat(
# 				(
# 					audio_obj.data.mfcc_coeffs,
# 					audio_obj.data.f0,
# 				)...)
# 		end,
# 		:age_set =>
# 			(x) -> begin
# 				get_full(x)
# 				return hcat(
# 					(
# 						# audio_obj.data.mel_spectrogram,
# 						# audio_obj.data.log_mel,
# 						audio_obj.data.mfcc_coeffs,
# 						audio_obj.data.mfcc_delta,
# 						audio_obj.data.mfcc_deltadelta,
# 						# audio_obj.data.spectral_centroid,
# 						# audio_obj.data.spectral_crest,
# 						# audio_obj.data.spectral_decrease,
# 						# audio_obj.data.spectral_entropy,
# 						# audio_obj.data.spectral_flatness,
# 						# audio_obj.data.spectral_flux,
# 						# audio_obj.data.spectral_kurtosis,
# 						# audio_obj.data.spectral_rolloff,
# 						# audio_obj.data.spectral_skewness,
# 						# audio_obj.data.spectral_slope,
# 						# audio_obj.data.spectral_spread,
# 						# audio_obj.data.f0,
# 					)...)
# 			end])

# 	if !isnothing(audio_obj)
# 		if haskey(func_call, feat)
# 			func_call[feat](audio_obj)
# 		else
# 			println("feature $feat not available.")
# 			return nothing
# 		end
# 	else
# 		return nothing
# 	end
# end

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

#------------------------------------------------------------------------------#
#                    delete when migration to dev is complete                  #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                              define structures                               #
#------------------------------------------------------------------------------#
function audio_setup(
	sr::Int64;

	# fft
	stft_length::Int64 = sr <= 8000 ? 256 : 512,
	win::AbstractVector{Float64} = Float64[],
	win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
	win_length::Int64 = stft_length, 					# standard setting: round(Int, 0.03 * sr)
	overlap_length::Int64 = round(Int, stft_length / 2), # standard setting: round(Int, 0.02 * sr)
	win_norm::Bool = false,
	spec_norm::Symbol = :power, #:power, :magnitude, :winpower, :winmagnitude

	# spectrum
	freq_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
	spectrum_type::Symbol = :power, 					# :power, :magnitude

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

		# fft
		stft_length = stft_length,
		win = win,
		win_type = win_type,
		win_length = win_length,
		overlap_length = overlap_length,
		win_norm = win_norm,
		spec_norm = spec_norm,

		# spectrum
		freq_range = freq_range,
		spectrum_type = spectrum_type,

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
	elseif length(x) < setup.stft_length
		@warn("length of sample is < than fft windows, skipped.")
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
function get_fft(audio_obj::AudioObj)
	if isempty(audio_obj.data.fft)
		get_fft!(audio_obj.setup, audio_obj.data)
	end
	return audio_obj.data.fft'
end

function get_lin_spec(audio_obj::AudioObj)
	if isempty(audio_obj.data.lin_spectrogram)
		if isempty(audio_obj.data.fft)
			get_fft(audio_obj)
		end
		lin_spectrogram!(audio_obj.setup, audio_obj.data)
	end
	return audio_obj.data.lin_spectrogram, audio_obj.data.lin_frequencies
end

function get_mel_spec(audio_obj::AudioObj)
	if isempty(audio_obj.data.lin_spectrogram)
		if isempty(audio_obj.data.fft)
			get_fft(audio_obj)
		end
		get_mel_spec!(audio_obj.setup, audio_obj.data)
	end
	return audio_obj.data.mel_spectrogram, audio_obj.data.mel_frequencies
end

function get_log_mel(audio_obj::AudioObj)
	if isempty(audio_obj.data.log_mel)
		if isempty(audio_obj.data.mel_spectrogram)
			if isempty(audio_obj.data.fft)
				get_fft(audio_obj)
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
			if isempty(audio_obj.data.fft)
				get_fft(audio_obj)
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
				if isempty(audio_obj.data.fft)
					get_fft(audio_obj)
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
	if audio_obj.setup.spectral_spectrum == :lin && isempty(audio_obj.data.lin_spectrogram)
		if isempty(audio_obj.data.fft)
			get_fft(audio_obj)
		end
		get_lin_spec(audio_obj)
	elseif audio_obj.setup.spectral_spectrum == :mel && isempty(audio_obj.data.mel_spectrogram)
		if isempty(audio_obj.data.fft)
			get_fft(audio_obj)
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
	get_fft(audio_obj)
	get_mel_spec(audio_obj)
	get_log_mel(audio_obj)
	get_mfcc(audio_obj)
	get_mfcc_delta(audio_obj)
	if audio_obj.setup.spectral_spectrum == :lin
		get_lin_spec(audio_obj)
	end
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
		:fft => get_fft,
		:lin => get_lin_spec,
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

function get_features(
	x::AbstractVector{T},
	sr::Int64,
	feat::Symbol = :full;
	kwargs...,
) where {T <: AbstractFloat}
	get_features(audio_obj(x, sr; kwargs...), feat)
end

function get_features(
	filepath::String,
	sr::Int64,
	feat::Symbol = :full;
	kwargs...,
)
	get_features(audio_obj(filepath, sr; kwargs...), feat)
end