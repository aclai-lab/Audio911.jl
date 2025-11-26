# ---------------------------------------------------------------------------- #
#                                AbstractInfo                                  #
# ---------------------------------------------------------------------------- #
abstract type AbstractInfo end

# ---------------------------------------------------------------------------- #
#                                AbstractFrame                                 #
# ---------------------------------------------------------------------------- #
"""
    AbstractFrame

Base type for windowing function implementations.
"""
abstract type AbstractFrame end

"""
    get_size(f::AbstractFrame) -> Int

Get the window size parameter from a `AbstractFrame` object.

# See also: [`AbstractFrame`](@ref), [`get_step`](@ref), [`get_overlap`](@ref)
"""
get_size(f::AbstractFrame) = error("get_size is not implemented for type $(typeof(f)).")

"""
    get_step(f::AbstractFrame) -> Int

Get the window step parameter from a `AbstractFrame` object.

# See also: [`AbstractFrame`](@ref), [`get_size`](@ref), [`get_overlap`](@ref)
"""
get_step(f::AbstractFrame) = error("get_step is not implemented for type $(typeof(f)).")

"""
    get_overlap(f::AbstractFrame) -> Int

Get the overlap length between consecutive windows from a `AbstractFrame` object.

# See also: [`AbstractFrame`](@ref), [`get_size`](@ref), [`get_step`](@ref)
"""
get_overlap(f::AbstractFrame) = error("get_overlap is not implemented for type $(typeof(f)).")

"""
    get_data(f::AbstractFrame) -> Vector{<:AudioFormat}

Extract the buffered audio frames from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref), [`get_winframes`](@ref), [`get_window`](@ref)
"""
get_data(f::AbstractFrame) = error("get_data is not implemented for type $(typeof(f)).")

"""
    get_window(f::AbstractFrame) -> Vector{Float64}

Extract the window from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref), [`get_data`](@ref), [`get_winframes`](@ref)
"""
get_window(f::AbstractFrame) = error("get_window is not implemented for type $(typeof(f)).")

"""
    get_winframes(f::AbstractFrame) -> Matrix

Get the windowed frames with the window function applied element-wise.

# See also: [`AbstractFrame`](@ref), [`get_data`](@ref), [`get_window`](@ref)
"""
get_winframes(f::AbstractFrame) = error("get_winframes is not implemented for type $(typeof(f)).")

"""
    get_winsize(f::AbstractFrame) -> Int

Get the length of the window function from an `AbstractFrame` container.

# See also: [`AudioFrames`](@ref), [`get_window`](@ref), [`get_winframes`](@ref)
"""
get_winsize(f::AbstractFrame) = error("get_winsize is not implemented for type $(typeof(f)).")

"""
    get_info(f::AbstractFrame) -> NamedTuple

Extract metadata from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref)
"""
get_info(f::AbstractFrame) = error("get_info is not implemented for type $(typeof(f)).")

# ---------------------------------------------------------------------------- #
#                             AbstractSpectrogram                              #
# ---------------------------------------------------------------------------- #
"""
    AbstractSpectrogram

Abstract supertype for all spectrogram representations in the Audio911.jl package.

This type serves as the base for concrete spectrogram implementations that store
time-frequency representations of audio signals. Subtypes should implement
methods for accessing spectral data, frequency information, and metadata.

# Concrete Subtypes
- [`Stft`](@ref): Short-Time Fourier Transform spectrogram representation

# Interface
Concrete subtypes are expected to provide:
- Spectral data matrix (typically complex or real-valued)
- Frequency vector corresponding to spectral bins
- Metadata about the analysis parameters

# See also
[`Stft`](@ref), [`get_stft`](@ref), [`get_freq`](@ref), [`get_info`](@ref)
"""
abstract type AbstractSpectrogram end

"""
    get_spec(s::AbstractSpectrogram) -> Matrix{T}

Extract the spectrogram matrix from an `AbstractSpectrogram` container.

# See also: [`get_freq`](@ref), [`get_info`](@ref)
"""
get_data(s::AbstractSpectrogram) = error("get_spec is not implemented for type $(typeof(s)).")

"""
    get_freq(s::AbstractSpectrogram) -> Vector{Float64}

Extract the frequency vector from an `AbstractSpectrogram` container.

# See also: [`get_spec`](@ref), [`get_info`](@ref)
"""
get_freq(s::AbstractSpectrogram) = error("get_freq is not implemented for type $(typeof(s)).")

"""
    get_info(s::AbstractSpectrogram) -> NamedTuple

Extract metadata from an `AbstractSpectrogram` container.

# See also: [`get_spec`](@ref), [`get_freq`](@ref)
"""
get_info(s::AbstractSpectrogram) = error("get_info is not implemented for type $(typeof(s)).")