module SoleAudio

using DSP
using FFTW
using LinearAlgebra
using Parameters
using SpecialFunctions
using Statistics

# structures
export signal_setup, signal_data
# main functions
export takeFFT, lin_spectrogram, mel_spectrogram, _mfcc, spectral_features, f0

include("signalDataStructure.jl")
# windowing
include("windowing/windows.jl")
include("windowing/windowing.jl")
# fft
include("fft/conv.jl")
include("fft/f0.jl")
include("fft/fft.jl")
include("fft/lin.jl")
include("fft/mel.jl")
include("fft/spectral.jl")

end # module SoleAudio
