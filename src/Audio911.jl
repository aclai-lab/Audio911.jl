module Audio911
using  Reexport

# ---------------------------------------------------------------------------- #
#                                audio reader                                  #
# ---------------------------------------------------------------------------- #
@reexport using AudioReader: @format_str, File, AudioFormat, AudioFile, load
@reexport using AudioReader: get_origin_sr, get_nchannels, is_norm
export get_data, get_sr
using  AudioReader: AudioFormat, _convert_mono
import AudioReader: get_data, get_sr

# ---------------------------------------------------------------------------- #
#                           audio related packages                             #
# ---------------------------------------------------------------------------- #
using  FFTW, DSP

# ---------------------------------------------------------------------------- #
#                              external packages                               #
# ---------------------------------------------------------------------------- #
using  DataTreatments
@reexport using DataTreatments: movingwindow

# ---------------------------------------------------------------------------- #
#                                   types                                      #
# ---------------------------------------------------------------------------- #
abstract type AbstractInfo end

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
get_low(fr::FreqRange)   = fr.low
get_hi(fr::FreqRange)    = fr.hi

export get_freqs, get_low, get_hi
export FreqRange

# ---------------------------------------------------------------------------- #
#                                  modules                                     #
# ---------------------------------------------------------------------------- #
# reexport DSP's window functions
@reexport using DSP: rect, hanning, hamming, cosine, lanczos, triang
@reexport using DSP: bartlett, bartlett_hann, blackman

export AbstractFrames
export AudioFrames
export get_size, get_step, get_overlap
export get_window, get_winframes, get_winsize
export nchannels
export get_info
include("frames.jl")

export AbstractSpectrogram
export get_spec, get_freq, get_info
include("spectrograms/types.jl")

export Stft
export power, magnitude
include("spectrograms/stft.jl")

export LinSpec
export none, winpower, winmagnitude
include("spectrograms/lin_spec.jl")

# export get_melspec
# include("mel_spec.jl")

end
