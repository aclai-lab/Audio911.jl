using DSP
using FFTW
using Parameters
using PyCall

include("../src/windowing/windowing.jl")
include("../src/windowing/windows.jl")
include("../src/signalDataStructure.jl")
include("../src/fft/fft.jl")
include("../src/fft/lin.jl")
include("../src/fft/mel.jl")
include("../src/fft/spectral.jl")
include("../src/fft/f0.jl")

af = pyimport("audioflux")
librosa = pyimport("librosa")

sr_src = 16000
x, sr = librosa.load("/home/riccardopasini/Documents/Aclai/Datasets/SpcDS/SpcDS_gender_1000_60_100/WavFiles/common_voice_en_23616312.wav", sr=sr_src, mono=true)
FFTLength = 256
mel_num = 26

# setup and data structures definition
setup = signal_setup(
    sr=sr,

    # fft
    window_type=[:hann, :periodic],
    window_length=FFTLength,
    overlap_length=round(Integer, FFTLength * 0.500),
    window_norm=:false,

    # spectrum
    frequency_range=[0, round(Integer, sr / 2)],
    spectrum_type=:power,

    # mel
    mel_style=:htk,
    num_bands=mel_num,
    filterbank_design_domain=:linear,
    filterbank_normalization=:bandwidth,
    frequency_scale=:mel,

    # mfcc
    num_coeffs=13,
    rectification=:log,
    log_energy_pos=:replace,
    delta_window_length=9,
    delta_axe=2,

    # spectral
    spectral_spectrum=:linear
)

# convert to Float64
x = Float64.(x)
# preemphasis
# zi = 2 * x[1] - x[2]
# filt!(x, [1.0, -0.97], 1.0, x, [zi])
# normalize
# x = x ./ maximum(abs.(x))

data = signal_data(
    x=x
)

takeFFT(data, setup)
# lin_spectrogram(data, setup)
mel_spectrogram(data, setup)
# _mfcc(data, setup)
# spectral_features(data, setup)
# f0(data, setup)