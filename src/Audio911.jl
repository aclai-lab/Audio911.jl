module Audio911

using FFTW, DSP
using LinearAlgebra
using SpecialFunctions
using StatsBase
using Statistics, Roots
using NaNStatistics
using Polynomials
using Plots

using PyCall

function __init__()
    py"""
    import librosa as librosa
    import soundfile as soundfile

    def load_audio(fname, sr):
        x, sr_def = librosa.load(fname, sr=sr, mono=True)
        return x, sr_def

    def save_audio(fname, x, sr):
        soundfile.write(fname, x, samplerate=sr, subtype='PCM_16')
    """
end

include("structs/audio.jl")
export Audio, load_audio

include("windowing/windowing.jl")
include("windowing/windows.jl")
include("structs/stft.jl")
export Stft, get_stft

include("structs/lin_spec.jl")
export LinSpec, get_linspec

include("structs/mel_fbank.jl")
include("structs/cwt_fbank.jl")
export MelFbank, get_melfb
export CwtFbank, get_cwtfb

include("structs/mel_spec.jl")
export MelSpec, get_melspec

include("structs/mfcc.jl")
export Mfcc, get_mfcc
export Deltas, get_deltas

include("structs/spectral.jl")
export Spectral, get_spectrals

include("structs/f0.jl")
export F0, get_f0

include("structs/cwt.jl")
export Cwt, get_cwt

include("utils/histogram.jl")
include("utils/speech_detector.jl")
export speech_detector

end