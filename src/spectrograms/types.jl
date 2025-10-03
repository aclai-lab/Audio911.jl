# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
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

#------------------------------------------------------------------------------#
#                            spectrogram methods                               #
#------------------------------------------------------------------------------#
"""
    get_spec(s::AbstractSpectrogram) -> Matrix{T}

Extract the spectrogram matrix from an `AbstractSpectrogram` container.

# See also: [`get_freq`](@ref), [`get_info`](@ref)
"""
get_spec(s::AbstractSpectrogram) = error("get_spec is not implemented for type $(typeof(s)).")

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