module Audio911

using FFTW, DSP
using LinearAlgebra
using Parameters
using SpecialFunctions
using StatsBase
using Statistics, Roots
using NaNStatistics
using Polynomials

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
# include("audioFeaturesExtractor.jl")
# windowing
include("windowing/cswindows.jl")
include("windowing/windows.jl")
include("windowing/windowing.jl")
# fft
include("fft/conv.jl")
# include("fft/f0.jl")
# include("fft/fft.jl")
# include("fft/lin.jl")
# include("fft/mel.jl")
# include("fft/spectral.jl")
include("fft/stft.jl")
include("fft/filterbanks.jl")
# utils
include("utils/histogram.jl")
include("utils/speech_detector.jl")
include("utils/in_out.jl")
include("utils/trimaudio.jl")
# wavelets
include("wavelet/wavelets_data_structures.jl")
include("wavelet/waveletData.jl")
include("wavelet/cwt.jl")
include("wavelet/expand.jl")
include("wavelet/wpdec.jl")
include("wavelet/wpspectrum.jl")
# constant-q transform
# include("cqt/cqt.jl")

# structures
export AudioObj, Audio, Stft, LinSpec, Fbank
# audio features
export audio_obj, get_features

# utility functions
export speech_detector
export load_audio, save_audio, trim_audio, normalize_audio
# wavelets
export cwt, cwt_windowing, wpspectrum
# get_stft!
# TODO patch
extractfeatures = 
export extractfeatures

end # module Audio911
