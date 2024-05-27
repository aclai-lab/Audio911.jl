
# -------------------------------------------------------------------------- #
#                                   debug                                    #
# -------------------------------------------------------------------------- #
using Revise, Plots
using Audio911, Parameters, FFTW, DSP, StatsBase, NaNStatistics
using Unitful, NamedArrays

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

sr_src = 16000
x, sr = load_audio(wavfile, sr_src)
x = Float64.(x)

include("/home/paso/.julia/dev/Audio911.jl/src/windowing/windows.jl")
include("/home/paso/.julia/dev/Audio911.jl/src/windowing/windowing.jl")
# extra parameters:
stft_length = round(Int, 0.05 * sr)
win_type = (:hann, :periodic)
win_length = round(Int, 0.03 * sr)
overlap_length = round(Int, 0.02 * sr)
freq_range = (100, 500)
spec_norm = :winpower
apply_log = false
thresh = -1000

stft_spec, stft_freq = get_stft(x, sr; stft_length=stft_length, win_type=win_type, win_length=win_length, overlap_length=overlap_length, freq_range=freq_range, spec_norm=spec_norm, apply_log=apply_log)