module Audio911
using  Reexport

# ---------------------------------------------------------------------------- #
#                                audio reader                                  #
# ---------------------------------------------------------------------------- #
@reexport using AudioReader: @format_str, File, AudioFormat, AudioFile, load
@reexport using AudioReader: get_origin_sr, get_nchannels, is_norm
# export get_data, get_sr
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
"""
    Maybe{T}

Type alias for `Union{T, Nothing}`.
"""
const Maybe{T} = Union{T, Nothing}

export AbstractInfo

export AbstractFrames
export get_data
export get_size, get_step, get_overlap
export get_window, get_winframes, get_winsize
export get_info

export AbstractSpectrogram
export get_freq, get_sr, get_nfft, get_spectype
export get_windows

export AbstractFBank
export get_bandwidth
export get_nbands, get_scale, get_norm
export get_freqrange, get_semitonerange
include("types.jl")

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
include("frames.jl")

export Stft
export power, magnitude
include("fft/stft.jl")

export area, bandwidth
export FBank
include("fft/fbank.jl")

export LinSpec
export winpower, winmagnitude
include("fft/lin_spec.jl")

# export get_melspec

end
