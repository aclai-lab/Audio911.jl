using DataFrames
using StatsBase, Random, MLJ
using Audio911

# ---------------------------------------------------------------------------- #
#                                load audiofile                                #
# ---------------------------------------------------------------------------- #
filename = "test.wav"
filepath = joinpath(@__DIR__, "test_files", filename)

audio = Audio911.load(filepath; norm=false, mono=true)
sr = samplerate(audio)

win = Audio911.MovingWindow(window_size=481, window_step=160)
frames = get_frames(audio; win, type=hamming)
frame1 = get_frames(audio; win, type=hamming, periodic=false)

stft = get_stft(frames; spectrum_type=:magnitude)
stft = get_stft(frames; frequency_range=(100,8000), spectrum_type=:magnitude)
mel_spec = get_melspec(stft)
