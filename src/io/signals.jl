# this code is taken from the package SampledSignals.jl
# (https://github.com/JuliaAudio/SampledSignals.jl)
# Copyright (c) 2015: Spencer Russell.

# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
"""
Represents a source of samples, such as an audio file, microphone input, or
SDR Receiver.

Subtypes should implement the `samplerate`, `nchannels`, `eltype`, and
`unsafe_read!` methods. `unsafe_read!` can assume that the samplerate, channel
count, and element type are all matching.
"""
abstract type AbstractSampleSource end
abstract type AbstractSampleBuf{T,N} <: AbstractArray{T,N} end

"""
Represents a multi-channel regularly-sampled buffer that stores its own sample
rate (in samples/second). The wrapped data is an N-dimensional array. A 1-channel
sample can be represented with a 1D array or an Mx1 matrix, and a C-channel
buffer will be an MxC matrix. So a 1-second stereo audio buffer sampled at
44100Hz with 32-bit floating-point samples in the time domain would have the
type SampleBuf{Float32, 2}.
"""
mutable struct SampleBuf{T,N} <: AbstractSampleBuf{T,N}
    data :: Array{T,N}
    sr   :: Int

    SampleBuf(data::Array{T,N}, sr::Int) where {T,N} = new{T,N}(data, sr)
end

