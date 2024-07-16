module Audio911

using FFTW, DSP
using LinearAlgebra
using SpecialFunctions
using StatsBase
using Statistics, Roots
using NaNStatistics
using Polynomials

using Parameters # DA CANCELLARE
using PyCall

function __init__()
    py"""
    import librosa as librosa
    import soundfile as soundfile

    def load_audio(filename, sr):
        x, sr_def = librosa.load(filename, sr=sr, mono=True)
        return x, sr_def

    def save_audio(filename, x, sr):
        soundfile.write(filename, x, samplerate=sr, subtype='PCM_16')
    """
end

include("signalDataStructure.jl")
include("spectrograms.jl")
# include("audioFeaturesExtractor.jl")
# windowing
include("windowing/cswindows.jl")
include("windowing/windows.jl")
include("windowing/windowing.jl")
# fft
include("fft/conv.jl")
# include("fft/f0.jl")
# include("fft/mel.jl")
# include("fft/spectral.jl")
include("fft/stft.jl")
# include("fft/lin.jl")
include("fft/fbank.jl")
# utils
include("utils/histogram.jl")
include("utils/speech_detector.jl")
include("utils/in_out.jl")
# include("utils/trimaudio.jl")
# wavelets
include("wavelets/wavelets_data_structures.jl")
include("wavelets/waveletData.jl")
include("wavelets/cwt_fbank.jl")
include("wavelets/cwt.jl")
include("wavelets/expand.jl")
include("wavelets/wpdec.jl")
include("wavelets/wpspectrum.jl")
# constant-q transform
# include("cqt/cqt.jl")

# structures
export AudioRack
export Audio, Stft, LinSpec, Fbank

# audio features
export get_spectrogram!
# stft based
export get_stft!, get_fbank!
# wavelets based
export get_cwt_fb!, get_cwt!

# utility functions
export speech_detector
export load_audio, save_audio, trim_audio, normalize_audio
# wavelets
# export cwt, cwt_windowing, wpspectrum
# get_stft!
# TODO patch
extractfeatures = 
export extractfeatures

end # module Audio911
