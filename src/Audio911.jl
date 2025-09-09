module Audio911
using  Reexport

# ---------------------------------------------------------------------------- #
#                                audio reader                                  #
# ---------------------------------------------------------------------------- #
@reexport using AudioReader: @format_str, File, AudioFormat, AudioFile, load
@reexport using AudioReader: data, sr
import AudioReader: nchannels

# ---------------------------------------------------------------------------- #
#                           audio related packages                             #
# ---------------------------------------------------------------------------- #
using FFTW

# ---------------------------------------------------------------------------- #
#                              external packages                               #
# ---------------------------------------------------------------------------- #
using  SoleBase: movingwindow

# ---------------------------------------------------------------------------- #
#                                   types                                      #
# ---------------------------------------------------------------------------- #
"""
    Maybe{T}

Type alias for `Union{T, Nothing}`.
"""
const Maybe{T} = Union{T, Nothing}

# ---------------------------------------------------------------------------- #
#                              frequency range                                 #
# ---------------------------------------------------------------------------- #
abstract type AbstractRange end

struct FreqRange{T} <: AbstractRange
    low :: T
    hi  :: T

    function FreqRange(low::T, hi:: T) where T
        new{T}(low, hi)
    end
end

Base.eltype(::Type{FreqRange{T}}) where T = T
Base.collect(fr::FreqRange) = collect(fr.low:fr.hi)

get_freqs(fr::FreqRange) = (fr.low, fr.hi)


export FreqRange

# ---------------------------------------------------------------------------- #
#                                  modules                                     #
# ---------------------------------------------------------------------------- #
export AbstractWinFunction, WinFunction, MovingWindow
include("windowing/windows.jl")
export AudioFrames
export get_frames, get_wsize, get_wstep, get_ovrlap, nchannels
include("windowing/audioframes.jl")

export Stft
export get_stft, get_stft_freq, get_info
include("fft/stft.jl")

export get_mel_spec
include("fft/mel_spec.jl")

# using FFTW, DSP
# using LinearAlgebra
# using Parameters
# using SpecialFunctions
# using StatsBase
# using Statistics, Roots
# using NaNStatistics

# include("signal_structure.jl")
# include("audioFeaturesExtractor.jl")
# # windowing
# include("windowing/cswindows.jl")
# include("windowing/windows.jl")

# # fft
# include("fft/conv.jl")
# include("fft/f0.jl")
# include("fft/lin.jl")

# include("fft/spectral.jl")
# # utils
# include("utils/histogram.jl")
# include("utils/speech_detector.jl")
# include("utils/in_out.jl")
# include("utils/trimaudio.jl")
# # wavelets
# include("wavelet/cwt.jl")
# # constant-q transform
# include("cqt/cqt.jl")

# # structures
# export AudioSetup, AudioData, AudioObj
# # audio features
# export audio_obj, get_features

# # utility functions
# export speech_detector
# export load_audio, save_audio, trim_audio, normalize_audio
# # wavelets
# export cwt, cwt_windowing
# get_fft!
# # TODO patch
# extractfeatures = 
# export extractfeatures

# # modular
# export get_frames, get_frames!, get_stft, get_stft!

# export get_frames2, get_stft2

end # module Audio911
