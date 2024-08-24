using Pkg
Pkg.activate("/home/paso/Documents/Aclai/audio-rules2024")
using Revise, Audio911, BenchmarkTools
using StaticArrays

using SpecialFunctions, Roots
using Statistics

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

# --- audio ------------------------------------------------------------------ #
sr = 16000
audio = load_audio(file=wavfile, sr=sr, norm=true);

# --- result ----------------------------------------------------------------- #
# result = audio911("stft(fft_length=1024), melfb(type=mel), cwt(source=mel_fb), mfcc(source=cwt, ncoeffs=16)")