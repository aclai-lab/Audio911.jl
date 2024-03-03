"""
    SoleAudio signal data structures

    uses package Parameter for @with_kw mutable struct

    signal_setup stores all datas that has to be shared in SoleAudio module
    signal_data stores all results from signal analysis
"""
@with_kw mutable struct signal_setup
    sr::Int64

    # fft
    fft_length::Int64 = 0
    window_type::Vector{Symbol} = [:hann, :periodic]
    window_length::Int64 = 0
    overlap_length::Int64 = 0
    window_norm::Bool = true

    # spectrum
    frequency_range::Vector{Int64} = []
    lin_frequencies::Vector{Float64} = []
    band_edges::AbstractVector{AbstractFloat} = []
    spectrum_type::Symbol = :power # :power, :magnitude

    # mel
    mel_style::Symbol = :htk # :htk, :slaney
    mel_bands::Int64 = 26
    mel_frequencies::Vector{Float64} = []
    filterbank_design_domain::Symbol = :linear
    filterbank_normalization::Symbol = :bandwidth # :bandwidth, :area, :none
    frequency_scale::Symbol = :mel # TODO :mel, :bark, :erb

    # mfcc
    num_coeffs::Int64 = 13
    normalization_type::Symbol = :standard # :standard, :dithered
    rectification::Symbol = :log
    log_energy_source::Symbol = :standard # :standard (after windowing), :mfcc
    log_energy_pos::Symbol = :append #:append, :replace, :none
    delta_window_length::Int64 = 9
    delta_matrix::Symbol = :standard # :standard, :transposed

    # spectral
    spectral_spectrum::Symbol = :linear # :linear, :mel
end

@with_kw mutable struct signal_data
    x::AbstractVector{Float64} = []

    # fft
    fft::AbstractArray{Float64} = []
    fft_window::Vector{Float64} = []

    # linear spectrum
    lin_spectrogram::AbstractArray{Float64} = []
    lin_frequencies::Vector{Float64} = []

    # mel_spectrum
    mel_filterbank::AbstractArray{Float64} = []
    mel_spectrogram::AbstractArray{Float64} = []

    # mfcc
    mfcc_coeffs::AbstractArray{Float64} = []
    mfcc_delta::AbstractArray{Float64} = []
    mfcc_deltadelta::AbstractArray{Float64} = []
    log_energy::Vector{Float64} = []

    # f0
    f0::Vector{Float64} = []

    # spectral
    spectral_centroid::Vector{Float64} = []
    spectral_crest::Vector{Float64} = []
    spectral_decrease::Vector{Float64} = []
    spectral_entropy::Vector{Float64} = []
    spectral_flatness::Vector{Float64} = []
    spectral_flux::Vector{Float64} = []
    spectral_kurtosis::Vector{Float64} = []
    spectral_rolloff::Vector{Float64} = []
    spectral_skewness::Vector{Float64} = []
    spectral_slope::Vector{Float64} = []
    spectral_spread::Vector{Float64} = []
end

# The other alternative, which is also common, is to define a structure with all parameters and pass that instead, like:

# julia> Base.@kwdef struct Options
#            a::Int = 1
#            b::Int = 2
#        end
# Options

# julia> _f(x; opt = Options()) = x + opt.a + opt.b
# _f (generic function with 1 method)

# julia> f(x; opt=Options()) = _f(x; opt=opt)
# f (generic function with 1 method)

# julia> f(1)
# 4

# julia> f(1; opt = Options(a=2))
# 5