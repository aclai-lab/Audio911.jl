using Audio911

using SpecialFunctions
using Statistics, Roots
using FFTW
using Parameters
using Plots
include("/home/paso/.julia/dev/Audio911.jl/src/wavelet/cwt.jl")

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")

sr_setup = 8000
x, sr = load_audio("$TESTPATH/common_voice_en_23616312.wav", sr=sr_setup)

window_length = 256
frequency_range=(80, 3000)
# mel_bands = 26
# num_coeffs = 13

setup = AudioSetup(
    sr=sr,
    # fft
    window_type=(:hann, :periodic),
    window_length=window_length,
    overlap_length=Int(round(window_length * 0.500)),
    window_norm=true,
    # spectrum
    frequency_range=frequency_range,
    spectrum_type=:power, # :power, :magnitude
    # mel
    mel_style=:htk, # :htk, :slaney
    mel_bands=mel_bands,
    filterbank_design_domain=:linear,
    filterbank_normalization=:bandwidth, # :bandwidth, :area, :none
    frequency_scale=:mel,
    # mfcc
    num_coeffs=num_coeffs,
    normalization_type=:dithered, # :standard, :dithered
    rectification=:log,
    log_energy_source=:standard, # :standard (after windowing), :mfcc
    log_energy_pos=:none, #:append, :replace, :none
    delta_window_length=9,
    delta_matrix=:standard, # :standard, :transposed
    # spectral
    spectral_spectrum=:mel # :linear, :linear_focused, :mel
)

data = AudioData(
    x=Float64.(x)
)

get_fft!(data, setup)

cwt_spectrum, _ = cwt(data.x, setup.sr, frequency_range=(80,3000))
cwt_spectrum = abs.(cwt_spectrum')

data.mel_spectrogram = cwt_windowing(cwt_spectrum, 32)
setup.mel_bands = setup.num_coeffs = size(data.mel_spectrogram, 2)

# mel_spectrogram(data, setup)
_mfcc(data, setup)
# lin_spectrogram(data, setup)
spectral_features(data, setup)
f0(data, setup)

# setup.frequency_range = (80, 1000)