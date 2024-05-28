"""
# optional arguments, default settings:
# stft_length::Int64 = 256,
# win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# win_length::Int64 = stft_length,
# overlap_length::Int64 = round(Int, stft_length / 2),
# window_norm::Bool = false,

# # spectrum
# freq_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
# spectrum_type::Symbol = :power,  # available options :power, :magnitude

# # mel
# mel_style::Symbol = :htk, # available options :htk, :slaney
# mel_bands::Int64 = 26,
# filterbank_design_domain::Symbol = :linear,
# filterbank_normalization::Symbol = :bandwidth, # available options :bandwidth, :area, :none
# frequency_scale::Symbol = :mel, # TODO :mel, :bark, :erb

# # mfcc
# mfcc_coeffs::Int64 = 13,
# normalization_type::Symbol = :dithered, # available options :standard, :dithered
# use_dct::Bool = true
# rectification::Symbol = :log, # available options :log, :cubic_root
# log_energy_source::Symbol = :standard, # available options :standard (after windowing), :mfcc
# log_energy_pos::Symbol = :none, # available options :append, :replace, :none
# delta_win_length::Int64 = 9,
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
using Plots
using Audio911

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

# -------------------------------------------------------------------------- #
#                                 parameters                                 #
# -------------------------------------------------------------------------- #
sr_src = 16000

# -------------------------------------------------------------------------- #
#                                 load audio                                 #
# -------------------------------------------------------------------------- #

x, sr = load_audio(wavfile, sr_src)

# -------------------------------------------------------------------------- #
#                        audio pre process utilities                         #
# -------------------------------------------------------------------------- #

# normalize audio signal
x = normalize_audio(x)
# speech detector trim
x, _ = speech_detector(x, sr)

# -------------------------------------------------------------------------- #
#     usage example 1: create audio object and call get features on that     #
# -------------------------------------------------------------------------- #
audio_1 = audio_obj(x, sr)

# full features
full_1 = get_features(audio_1, :full)

# single features
fft_1 = get_features(audio_1, :fft)
(lin_1, linfreq_1) = get_features(audio_1, :lin)
(mel_1, melfreq_1) = get_features(audio_1, :mel)
logmel_1 = get_features(audio_1, :logmel)
mfcc_1 = get_features(audio_1, :mfcc)
mfccdelta_1 = get_features(audio_1, :mfcc_delta)
spectral_1 = get_features(audio_1, :spectral)
f0_1 = get_features(audio_1, :f0)
cqt_1 = get_features(audio_1, :cqt)

# submit parameters
custom_1 = audio_obj(
	x,
	sr,
	stft_length = 1024,
	spectral_spectrum = :mel,
	freq_range = (50, 1000),
	mel_style = :slaney,
)
# -------------------------------------------------------------------------- #
#            usage example 2: calling features on preloaded audio            #
# -------------------------------------------------------------------------- #

# full features
full_2 = get_features(x, sr, :full)

# single features
fft_2 = get_features(x, sr, :fft)
(lin_2, linfreq_2) = get_features(x, sr, :lin)
(mel_2, melfreq_2) = get_features(x, sr, :mel)
logmel_2 = get_features(x, sr, :logmel)
mfcc_2 = get_features(x, sr, :mfcc)
mfccdelta_2 = get_features(x, sr, :mfcc_delta)
spectral_2 = get_features(x, sr, :spectral)
f0_2 = get_features(x, sr, :f0)
cqt_2 = get_features(x, sr, :cqt)

# submit parameters
custom_2 = get_features(
	x,
	sr,
	:full,
	stft_length = 1024,
	spectral_spectrum = :mel,
	freq_range = (50, 1000),
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

(mel_3, melfreq_3) = get_features(x, sr, :mel, frequency_scale = :chroma)

# -------------------------------------------------------------------------- #
#                        	   Gio's semitones                               #
# -------------------------------------------------------------------------- #
(tune_1, tunefreq_1) = get_features(x, sr, :mel, mel_style = :tuned, st_peak_range = (50, 500))

# -------------------------------------------------------------------------- #
#                           modular implementation                           #
# -------------------------------------------------------------------------- #

# --------------------------- stft standalone ------------------------------ #
# example 1:
stft_1, freq_1 = get_stft(x, sr)
heatmap(
	collect(1:size(stft_1, 2)) / sr, 
	freq_1, 
	stft_1, 
	xlabel="Time", 
	ylabel="Frequency", 
	title="Stft Example 1", color=:plasma, clims=(0, 20)
)

# example 2.1, with args, power spectrum
stft_2_1, freq_2_1 = get_stft(
	x,
	sr;
	stft_length = sr <= 8000 ? 256 : 512,
	win_type = (:hann, :periodic),
	win_length = sr <= 8000 ? 256 : 512,
	spec_norm = :power,
	freq_range = (0, floor(Int, sr / 2))
)
heatmap(
	collect(1:size(stft_2_1, 2)) / sr, 
	freq_2_1, 
	stft_2_1, 
	xlabel="Time", 
	ylabel="Frequency", 
	title="Stft Example 1", color=:plasma, clims=(0, 20)
)

# example 2.2, with args, magnitude spectrum
stft_2_2, freq_2_2 = get_stft(
	x,
	sr;
	stft_length = sr <= 8000 ? 256 : 512,
	win_type = (:hann, :periodic),
	win_length = sr <= 8000 ? 256 : 512,
	spec_norm = :magnitude,
	freq_range = (0, floor(Int, sr / 2))
)
heatmap(
	collect(1:size(stft_2_2, 2)) / sr, 
	freq_2_2, 
	stft_2_2, 
	xlabel="Time", 
	ylabel="Frequency", 
	title="Stft Example 1", color=:plasma, clims=(0, 20)
)

# example 2.3, with args, power spectrum, windows normalized
stft_2_3, freq_2_3 = get_stft(
	x,
	sr;
	stft_length = sr <= 8000 ? 256 : 512,
	win_type = (:hann, :periodic),
	win_length = sr <= 8000 ? 256 : 512,
	spec_norm = :winpower,
	freq_range = (0, floor(Int, sr / 2))
)
heatmap(
	collect(1:size(stft_2_3, 2)) / sr, 
	freq_2_3, 
	stft_2_3, 
	xlabel="Time", 
	ylabel="Frequency", 
	title="Stft Example 1", color=:plasma, clims=(0, 20)
)

# example 2.4, with args, magnitude spectrum, windows normalized
stft_2_4, freq_2_4 = get_stft(
	x,
	sr;
	stft_length = sr <= 8000 ? 256 : 512,
	win_type = (:hann, :periodic),
	win_length = sr <= 8000 ? 256 : 512,
	spec_norm = :winmagnitude,
	freq_range = (0, floor(Int, sr / 2))
)
heatmap(
	collect(1:size(stft_2_4, 2)) / sr, 
	freq_2_4, 
	stft_2_4, 
	xlabel="Time", 
	ylabel="Frequency", 
	title="Stft Example 1", color=:plasma, clims=(0, 20)
)

# example 3, with stft larger than window frames length, for a better frequency resolution
stft_3, freq_3 = get_stft(
	x,
	sr;
	stft_length = sr <= 8000 ? 2048 : 4096,
	win_length = sr <= 8000 ? 256 : 512,
	spec_norm = :power,
)
heatmap(
	collect(1:size(stft_3, 2)) / sr, 
	freq_3, 
	stft_3, 
	xlabel="Time", 
	ylabel="Frequency", 
	title="Stft Example 1", color=:plasma, clims=(0, 20)
)

# --------------------------- stft audio object ------------------------------ #
# example 4, utilizing audio_obj
audio_1 = audio_objdev(x, sr)
get_stft!(audio_1)

# # example 5, utilizing audio_obj with args
# audio_5 = audio_obj(x,
# 	sr,
# 	stft_length = 512,
# 	spec_norm = :magnitude,
# 	win_type = (:hamming, :symmetric),
# 	freq_range = (100, 3000))
# get_stft!(audio_5)

# # example 6, utilizing audio_obj modular
# audio_6 = audio_obj(x, sr)
# get_frames!(audio_6)
# get_stft!(audio_6)
# display(audio_6.data.stft)

# ################
# framest = get_frames2(x)
# stft1, stft_freq1 = get_stft2(framest, sr)
# fft22 = get_features(x, sr, :fft)