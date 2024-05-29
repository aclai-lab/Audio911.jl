using Revise
using Plots
using Audio911
using StatsBase

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
#          usage: create audio object and call get features on that          #
# -------------------------------------------------------------------------- #
audio_1 = audio_obj(x, sr)
# audio_1 = audio_obj(wavfile, sr)

# full features
full_1 = get_features(audio_1, :full)

# single features
fft_1 = get_features(audio_1, :stft)
mel_1, melfreq_1 = get_features(audio_1, :mel)
# logmel_1 = get_features(audio_1, :logmel)
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
############################# DEPRECATED #####################################
# full features
# full_2 = get_features(x, sr, :full)

# # single features
# fft_2 = get_features(x, sr, :fft)
# # (lin_2, linfreq_2) = get_features(x, sr, :lin)
# (mel_2, melfreq_2) = get_features(x, sr, :mel)
# logmel_2 = get_features(x, sr, :logmel)
# mfcc_2 = get_features(x, sr, :mfcc)
# mfccdelta_2 = get_features(x, sr, :mfcc_delta)
# spectral_2 = get_features(x, sr, :spectral)
# f0_2 = get_features(x, sr, :f0)
# cqt_2 = get_features(x, sr, :cqt)

# submit parameters
# custom_2 = get_features(
# 	x,
# 	sr,
# 	:full,
# 	stft_length = 1024,
# 	spectral_spectrum = :mel,
# 	freq_range = (50, 1000),
# 	mel_style = :slaney,
# 	mel_bands = 40,
# )

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
#                           modular implementation                           #
# -------------------------------------------------------------------------- #



# -------------------------------------------------------------------------- #
#                              features scaling                              #
# -------------------------------------------------------------------------- #
# features scaling references
# https://medium.com/@shivanipickl/what-is-feature-scaling-and-why-does-machine-learning-need-it-104eedebb1c9
# https://scikit-learn.org/stable/auto_examples/preprocessing/plot_scaling_importance.html

X = full_1

dt = fit(ZScoreTransform, X; dims=1)
zX = StatsBase.transform(dt, X)

dt = fit(UnitRangeTransform, X, dims=1)
uX = StatsBase.transform(dt, X)