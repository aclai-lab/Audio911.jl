"""
# optional arguments, default settings:
# fft_length::Int64 = 256,
# window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# window_length::Int64 = fft_length,
# overlap_length::Int64 = round(Int, fft_length / 2),
# window_norm::Bool = false,

# # spectrum
# frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
# spectrum_type::Symbol = :power,  # available options :power, :magnitude

# # mel
# mel_style::Symbol = :htk, # available options :htk, :slaney
# mel_bands::Int64 = 26,
# filterbank_design_domain::Symbol = :linear,
# filterbank_normalization::Symbol = :bandwidth, # available options :bandwidth, :area, :none
# frequency_scale::Symbol = :mel, # TODO :mel, :bark, :erb

# # mfcc
# num_coeffs::Int64 = 13,
# normalization_type::Symbol = :dithered, # available options :standard, :dithered
# use_dct::Bool = true
# rectification::Symbol = :log, # available options :log, :cubic_root
# log_energy_source::Symbol = :standard, # available options :standard (after windowing), :mfcc
# log_energy_pos::Symbol = :none, # available options :append, :replace, :none
# delta_window_length::Int64 = 9,
# delta_matrix::Symbol = :transposed, # available options :standard, :transposed

# # spectral
# spectral_spectrum::Symbol = :lin, # available options :lin, :mel

# # f0
# f0_method::Symbol = :nfc,
# f0_range::Tuple{Int64, Int64} = (50, 400),
# median_filter_length::Int64 = 1,

# constant-q transform	
# bins_octave::Int64 = 12, default is 12: semitone scale
# freq_limits::Tuple{Float64, Float64} = (sr / x_length, sr / 2),
# transform_type::Symbol = :full, # available options :full, :reduced, :sparse
"""

using Revise
using Audio911

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

# -------------------------------------------------------------------------- #
#                                 parameters                                 #
# -------------------------------------------------------------------------- #
sr_src = 8000

# -------------------------------------------------------------------------- #
#                                 load audio                                 #
# -------------------------------------------------------------------------- #

x, sr = load_audio(wavfile, sr = sr_src)

# -------------------------------------------------------------------------- #
#                        audio pre process utilities                         #
# -------------------------------------------------------------------------- #

# normalize audio signal
x = normalize_audio(x)
# speech detector trim
# x = speech_detector(x, sr) BROKEN!

# -------------------------------------------------------------------------- #
#     usage example 1: create audio object and call get features on that     #
# -------------------------------------------------------------------------- #
audio_1 = audio_obj(x, sr)

# full features
full_2 = get_features(audio_1, :full)

# single features
fft_2 = get_features(audio_1, :fft)
lin_2 = get_features(audio_1, :lin)
mel_2 = get_features(audio_1, :mel)
logmel_2 = get_features(audio_1, :logmel)
mfcc_2 = get_features(audio_1, :mfcc)
mfccdelta__2 = get_features(audio_1, :mfcc_delta)
spectral_2 = get_features(audio_1, :spectral)
f0_2 = get_features(audio_1, :f0)
cqt_2 = get_features(audio_1, :cqt)

# submit parameters
custom_1 = audio_obj(
	x,
	sr,
	fft_length = 1024,
	spectral_spectrum = :mel,
	frequency_range = (50, 1000),
	mel_style = :slaney,
)
# -------------------------------------------------------------------------- #
#            usage example 2: calling features on preloaded audio            #
# -------------------------------------------------------------------------- #

# full features
full_2 = get_features(x, sr, :full)

# single features
fft_2 = get_features(x, sr, :fft)
lin_2 = get_features(x, sr, :lin)
mel_2 = get_features(x, sr, :mel)
logmel_2 = get_features(x, sr, :logmel)
mfcc_2 = get_features(x, sr, :mfcc)
mfccdelta__2 = get_features(x, sr, :mfcc_delta)
spectral_2 = get_features(x, sr, :spectral)
f0_2 = get_features(x, sr, :f0)
cqt_2 = get_features(x, sr, :cqt)

# submit parameters
custom_2 = get_features(
	x,
	sr,
	:full,
	fft_length = 1024,
	spectral_spectrum = :mel,
	frequency_range = (50, 1000),
	mel_style = :slaney,
	mel_bands = 40,
)

# -------------------------------------------------------------------------- #
#                      usage example 3: submit filepath                      #
# -------------------------------------------------------------------------- #

sr = 8000
audio_3 = audio_obj(wavfile, sr)
audio_3_1 = get_features(wavfile, sr, :full)

# -------------------------------------------------------------------------- #
#                                 save audio                                 #
# -------------------------------------------------------------------------- #



# -------------------------------------------------------------------------- #
#              exceptions: audio sample is empty or too small                #
# -------------------------------------------------------------------------- #

bad_x = Float64[]
bad_1 = get_features(bad_x, sr)

bad_x = (rand(1.0:3.0, 128))
bad_1 = get_features(bad_x, sr)

mel_1 = get_features(x, sr, :mel, frequency_scale = :chroma)

# -------------------------------------------------------------------------- #
#                        	   Gio's semitones                               #
# -------------------------------------------------------------------------- #
semi_1 = get_features(x, sr, :mel, mel_style=:semitones)
