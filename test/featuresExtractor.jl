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
    # window_length::Int=Int(round(0.03 * sr)),
    # overlap_length::Int=Int(round(0.02 * sr)),
    window_norm=:false,

    # spectrum
    frequency_range=[0, round(Integer, sr / 2)],
    spectrum_type=:power,

    # mel
    mel_style=:htk,
    mel_bands=mel_num,
    filterbank_design_domain=:linear,
    filterbank_normalization=:bandwidth,
    frequency_scale=:mel,

    # mfcc
    num_coeffs=13,
    rectification=:log,
    # log_energy_source=:standard,
    log_energy_pos=:append,
    delta_window_length=9,
    delta_matrix=:transposed,

    # spectral
    spectral_spectrum=:linear
)

# convert to Float64
x = Float64.(x)

# preemphasis
# not siutable for our kind of experiments, maybe for speak recognition: needs to look over it.
# zi = 2 * audio[1] - audio[2]
# filt!(audio, [1.0, -0.97], 1.0, audio, [zi])
# normalize
# TODO
# already normalized. think to remove previous normalization and move it here.
# audio = audio ./ maximum(abs.(audio))

data = signal_data(
    x=x
)

takeFFT(data, setup)
# lin_spectrogram(data, setup)
mel_spectrogram(data, setup)
# _mfcc(data, setup)
# spectral_features(data, setup)
# f0(data, setup)

audio_features = vcat((
    data.mfcc_coeffs',
    data.mfcc_delta',
    # data.mfcc_deltadelta',
    data.mel_spectrogram[:, 1:13]'
    # data.mel_spectrogram'
)...)