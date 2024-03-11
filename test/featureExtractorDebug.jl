using DSP
using FFTW
using Parameters
using PyCall
using LinearAlgebra

# using SoleAudio

include("../src/windowing/windowing.jl")
include("../src/windowing/windows.jl")
include("../src/signalDataStructure.jl")
include("../src/fft/fft.jl")
include("../src/fft/lin.jl")
include("../src/fft/mel.jl")
include("../src/fft/spectral.jl")
include("../src/fft/f0.jl")
include("../src/utils/in_out.jl")


# af = pyimport_conda("audioflux")
librosa = pyimport("librosa")

sr_src = 8000
x, sr = librosa.load("/home/riccardopasini/.julia/dev/SoleAudio.jl/test/common_voice_en_23616312.wav", sr=sr_src, mono=true)
fft_length = 256
frequency_range=Int[0, sr/2]
mel_bands = 26
num_coeffs = 13

setup = signal_setup(
    sr=sr,
    # fft
    window_type=[:hann, :periodic],
    window_length=fft_length,
    overlap_length=Int(round(fft_length * 0.500)),
    window_norm=false,
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
    spectral_spectrum=:linear # :linear, :linear_focused, :mel
)

# convert to Float64
x = Float64.(x)

data = signal_data(
    x=x
)

takeFFT(data, setup)
mel_spectrogram(data, setup)
_mfcc(data, setup)
lin_spectrogram(data, setup)
spectral_features(data, setup)
f0(data, setup) # pay attention to fft length

data.spectral_entropy