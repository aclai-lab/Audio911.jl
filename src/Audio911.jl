module Audio911
using  Reexport

# ---------------------------------------------------------------------------- #
#                                audio reader                                  #
# ---------------------------------------------------------------------------- #
@reexport using AudioReader: @format_str, File, AudioFormat, AudioFile, load
@reexport using AudioReader: get_origin_sr, get_nchannels, is_norm

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
using Plots

# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
abstract type AbstractSetup       end
abstract type AbstractFrame       end
abstract type AbstractSpectrogram end
abstract type AbstractFBank       end

# ---------------------------------------------------------------------------- #
#                                   types                                      #
# ---------------------------------------------------------------------------- #
# type alias for `Union{T, Nothing}`
const Maybe{T} = Union{T, Nothing}

"""
    AudioData

Type alias for audio sample data. Represents a value that can be either `Float64` or `Float32`.

This alias is used for audio arrays and sample values throughout the package to ensure type stability.

# Example

```julia
x::Vector{AudioData} = rand(Float32, 1024)
y::AudioData = 0.5
```
"""
const  AudioData = Union{Float64, Float32}
export AudioData

"""
    FreqRange

Type alias for a frequency range, represented as a tuple `(min, max)` of integer values.

This alias is used to specify the minimum and maximum frequency (in Hz) for audio processing routines.

# Example

```julia
fr::FreqRange = (20, 20000)  # 20 Hz to 20 kHz
```
"""
const  FreqRange = Tuple{T, T} where {T<:Int64}
 
get_low(r::FreqRange) = r[1]
get_hi(r::FreqRange)  = r[2]
export FreqRange, get_low, get_hi

"""
    ScaleRange

Type alias for a perceptual scale range, represented as a tuple `(min, max)` of `AudioData` values.

This alias is used to specify the minimum and maximum values in perceptual scales (mel, bark, ERB, semitones)
for audio processing routines.

# Example

```julia
mel_range::ScaleRange = (0.0, 45.0)     # Mel scale range
erb_range::ScaleRange = (0.0f0, 43.0f0) # ERB scale range (Float32)
```
"""
const  ScaleRange  = Tuple{T, T} where {T<:AudioData}

get_low(r::ScaleRange) = r[1]
get_hi(r::ScaleRange)  = r[2]
export ScaleRange

# ---------------------------------------------------------------------------- #
#                           spectrum normalizations                            #
# ---------------------------------------------------------------------------- #
winpower(f, w)     = f / sum(w).^2
winmagnitude(f, w) = f / sum(w)

export winpower, winmagnitude

# ---------------------------------------------------------------------------- #
#                                  modules                                     #
# ---------------------------------------------------------------------------- #
# reexport DSP's window functions
@reexport using DSP: rect, hanning, hamming, cosine, lanczos, triang
@reexport using DSP: bartlett, bartlett_hann, blackman

export AbstractFrames
export Frames
include("frames.jl")

export Stft
export power, magnitude
include("fft/stft.jl")

export FBank
export htk, slaney, bark
export area, bandwidth, none_norm
export auditory_fbank, gammatone_fbank
include("fft/fbank.jl")

export LinSpec
include("fft/lin_spec.jl")

export MelSpec
include("fft/mel_spec.jl")

export Mfcc
export mlog, cubic_root
include("fft/mfcc.jl")

# ---------------------------------------------------------------------------- #
#                                  methods                                     #
# ---------------------------------------------------------------------------- #
# general
export get_data, get_setup

# stft related
export get_size, get_step, get_overlap
export get_window, get_winframes, get_winsize

# spectrogram related
export get_freq, get_sr, get_nfft, get_spectrum
export get_windows

# filterbank related
export get_bandwidth
export get_nbands, get_scale, get_norm
export get_freqrange, get_semitonerange

# ---------------------------------------------------------------------------- #
#                                   plots                                      #
# ---------------------------------------------------------------------------- #
export plot
include("plots.jl")

end
