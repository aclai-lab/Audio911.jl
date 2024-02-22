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

# af = pyimport("audioflux")
librosa = pyimport("librosa")

sr_src = 8000
x, sr = librosa.load("/home/riccardopasini/Documents/Aclai/Datasets/0c40e715_nohash_0.wav", sr=sr_src, mono=true)
FFTLength = 256
mel_num = 26

setup = signal_setup(
    sr=sr,

    # fft
    window_type=[:hann, :periodic],
    window_length=FFTLength,
    overlap_length=Int(round(FFTLength * 0.500)),
    # window_length=Int(round(0.03 * sr)),
    # overlap_length=Int(round(0.02 * sr)),
    window_norm=:false,

    # spectrum
    frequency_range=Int[0, sr/2],
    spectrum_type=:power,

    # mel
    mel_style=:htk,
    mel_bands=mel_num,
    filterbank_design_domain=:linear,
    filterbank_normalization=:bandwidth,
    frequency_scale=:mel,

    # mfcc
    num_coeffs=13,
    normalization_type=:dithered,
    rectification=:log,
    log_energy_source=:standard,
    log_energy_pos=:replace,
    delta_window_length=9,
    delta_matrix=:transposed,

    # spectral
    spectral_spectrum=:linear
)

# convert to Float64
x = Float64.(x)

data = signal_data(
    x=x
)

takeFFT(data, setup)
lin_spectrogram(data, setup)
# mel_spectrogram(data, setup)
# _mfcc(data, setup)
# spectral_features(data, setup)
f0(data, setup)