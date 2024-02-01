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
    frequency_range::Vector{Float64} = []
    lin_frequencies::Vector{Float64} = []
    band_edges::AbstractVector{AbstractFloat} = []
    spectrum_type::Symbol = :power # :power, :magnitude

    # mel
    mel_style::Symbol = :htk # :htk, :slaney
    mel_bands::Int64 = 32
    mel_frequencies::Vector{Float64} = []
    filterbank_design_domain::Symbol = :linear
    filterbank_normalization::Symbol = :bandwidth # :bandwidth, :area, :none
    frequency_scale::Symbol = :mel # :mel, :erb

    # mfcc
    num_coeffs::Int64 = 13
    dct_type::Int64 = 2 # from 1 to 4
    rectification::Symbol = :log
    # log_energy_source::Symbol = :standard # :standard (after windowing), :mfcc
    log_energy_pos::Symbol = :append #:append, :replace, :none
    delta_window_length::Int64 = 9
    delta_matrix::Symbol = :standard # :standard, :transposed

    # spectral
    spectral_spectrum::Symbol = :linear
end

@with_kw mutable struct signal_data
    x::AbstractArray{Float64} = []

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