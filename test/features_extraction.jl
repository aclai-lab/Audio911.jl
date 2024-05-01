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
"""

using Revise
using Audio911

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"

# -------------------------------------------------------------------------- #
#                                 parameters                                 #
# -------------------------------------------------------------------------- #
sr_src = 8000

# -------------------------------------------------------------------------- #
#                                 load audio                                 #
# -------------------------------------------------------------------------- #

x, sr = load_audio(joinpath(TESTPATH, TESTFILE), sr = sr_src)

# -------------------------------------------------------------------------- #
#                        audio pre process utilities                         #
# -------------------------------------------------------------------------- #

# normalize audio signal
x = normalize_audio(x)
# speech detector trim
# x = speech_detector(x, sr) BROKEN!

# -------------------------------------------------------------------------- #
#                        usage example 1: julia style                        #
# -------------------------------------------------------------------------- #

# full features
full_1 = get_features(x, sr)
full_1 = get_features(x, sr, :full)

# single features
fft_1 = get_features(x, sr, :fft)
lin_1 = get_features(x, sr, :lin)
mel_1 = get_features(x, sr, :mel)
logmel_1 = get_features(x, sr, :logmel)
mfcc_1 = get_features(x, sr, :mfcc)
spectral_1 = get_features(x, sr, :spectral)
f0_1 = get_features(x, sr, :f0)

# submit parameters
custom_1 = get_features(x, sr, fft_length = 1024, spectral_spectrum = :mel, frequency_range = (50, 1000), mel_style = :slaney)

# -------------------------------------------------------------------------- #
#                       usage example 2: object style                        #
# -------------------------------------------------------------------------- #

audio = audio_features_obj(x, sr, f0_range = (50,700))

audio.get_fft()
audio.get_lin_spec()
audio.get_mel_spec()
audio.get_log_mel()
audio.get_mfcc()
audio.get_spectrals()
audio.get_f0()
audio.get_features(:full)

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

bad_x = Float64[]
bad_1 = audio_features_obj(bad_x, sr)

bad_x = (rand(1.0:3.0, 128))
bad_1 = audio_features_obj(bad_x, sr)
