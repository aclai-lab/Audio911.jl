module Audio911

import Base: show
import Plots: plot

using FFTW, DSP
using Base.Threads
using LinearAlgebra
using SpecialFunctions
using StatsBase
using Statistics, Roots
using NaNStatistics
using Polynomials
using Plots
using SparseArrays, StaticArrays
using SplitApplyCombine

using PyCall

function __init__()
    py"""
    import librosa as librosa
    import soundfile as soundfile

    def load_audio(file, sr):
        x, sr_def = librosa.load(file, sr=sr)
        return x, sr_def

    def save_audio(file, x, sr):
        soundfile.write(file, x, samplerate=sr, subtype='PCM_16')
    """
end

include("structs/audio.jl")
export Audio, load_audio, save_audio

# include("windowing/windows.jl")
include("structs/stft.jl")
export Stft, get_stft

include("structs/lin_spec.jl")
export LinSpec, get_linspec

include("structs/mel_fbank.jl")
include("structs/cwt_fbank.jl")
export MelFb, get_melfb
export CwtFb, get_cwtfb

include("structs/cwt.jl")
export Cwt, get_cwt

include("structs/mel_spec.jl")
export MelSpec, get_melspec

include("structs/mfcc.jl")
include("structs/deltas.jl")
export Mfcc, get_mfcc
export Deltas, get_deltas

include("structs/spectral.jl")
export Spectral, get_spectrals

include("structs/f0.jl")
export F0, get_f0

include("utils/histogram.jl")
include("utils/speech_detector.jl")
export speech_detector

include("interface.jl")
export plot
export audio_features

# DEBUG!
export _get_frames

end