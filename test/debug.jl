
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

include("/home/paso/.julia/dev/Audio911.jl/src/signalDataStructure.jl")
include("/home/paso/.julia/dev/Audio911.jl/src/windowing/windows.jl")
include("/home/paso/.julia/dev/Audio911.jl/src/windowing/windowing.jl")
include("/home/paso/.julia/dev/Audio911.jl/src/fft/stft.jl")
# extra parameters:
stft_length = sr <= 8000 ? 256 : 512
win_type = (:hann, :periodic)
win_length = stft_length
overlap_length = round(Int, stft_length / 2)
freq_range = (300, floor(Int, sr / 2))
spec_norm = :power
apply_log = false
thresh = -1000

stft_spec, stft_freq = get_stft(x, sr; stft_length=stft_length, win_type=win_type, win_length=win_length, overlap_length=overlap_length, freq_range=freq_range, spec_norm=spec_norm, apply_log=apply_log)