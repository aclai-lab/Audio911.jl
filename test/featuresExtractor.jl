using DSP

include("../src/signalDataStructure.jl")
include("../src/fft/fft.jl")
include("../src/fft/lin.jl")
include("../src/fft/mel.jl")
include("../src/fft/spectral.jl")
include("../src/fft/f0.jl")

function audioFeaturesExtraction(
    x::AbstractArray{T},
    sr::Int64;
    # default options
    # fft
    window_type::Symbol=:hann,
    window_length::Int=Int(round(0.03 * sr)),
    overlap_length::Int=Int(round(0.02 * sr)),

    # spectrum
    frequency_range::Vector{Int64}=[0, Int(round(sr / 2))],
    spectrum_type::Symbol=:power,

    # mel
    mel_style::Symbol=:default_htk, # :htk, :slaney, :default_htk, :default_slaney
    num_bands::Int=32,
    filterbank_design_domain::Symbol=:linear,
    filterbank_normalization::Symbol=:bandwidth,

    # mfcc
    num_coeffs::Int=13,
    rectification::Symbol=:log,
    log_energy_pos::Symbol=:append,
    delta_window_length::Int=9,

    # spectral
    spectral_spectrum::Symbol=:linear
) where {T<:AbstractFloat}

    # set default values for default :htk or :slaney, assuming sr = 16000
    # cita la fonte... comparative mfcc
    if (mel_style == :default_htk)
        num_bands = 24
        frequency_range = [0, 8000]
    elseif (mel_style == :default_slaney)
        num_bands = 40
        frequency_range = [133, 6854]
    end

    # setup and data structures definition
    setup = signal_setup(
        sr=sr,

        # fft
        window_type=window_type,
        window_length=window_length,
        overlap_length=overlap_length,

        # spectrum
        frequency_range=frequency_range,
        spectrum_type=spectrum_type,

        # mel
        mel_style=mel_style,
        num_bands=num_bands,
        filterbank_design_domain=filterbank_design_domain,
        filterbank_normalization=filterbank_normalization,

        # mfcc
        num_coeffs=num_coeffs,
        rectification=rectification,
        log_energy_pos=log_energy_pos,
        delta_window_length=delta_window_length,

        # spectral
        spectral_spectrum=spectral_spectrum
    )

    # convert to Float64
    x = Float64.(x)
    # preemphasis
    zi = 2 * x[1] - x[2]
    filt!(x, [1.0, -0.97], 1.0, x, [zi])
    # normalize
    x = x ./ maximum(abs.(x))
    # avoid 0.0 values
    # x[x.==0.0] .+= eps(0.0)

    data = signal_data(
        x=x
    )

    takeFFT(data, setup)
    lin_spectrogram(data, setup)
    mel_spectrogram(data, setup)
    mfcc(data, setup)
    spectral_features(data, setup)

    vcat(
        data.mel_spectrogram',
        data.mfcc_coeffs',
        data.mfcc_delta',
        data.mfcc_deltadelta',
        data.spectral_centroid',
        data.spectral_crest',
        data.spectral_decrease',
        data.spectral_entropy',
        data.spectral_flatness',
        data.spectral_flux',
        data.spectral_kurtosis',
        data.spectral_rolloff',
        data.spectral_skewness',
        data.spectral_slope',
        data.spectral_spread'
    )
end

# # debug
# using PyCall
# librosa = pyimport("librosa")
# sr_src = 16000
# x, sr = librosa.load("/home/riccardopasini/Documents/Aclai/Datasets/SpcDS/SpcDS_gender_1000_60_100/WavFiles/common_voice_en_23616312.wav", sr=sr_src, mono=true)

# # fft
# window_type = :hann
# window_length = Int(round(0.03 * sr))
# overlap_length = Int(round(0.02 * sr))
# # window_length = 256 # audioflux
# # overlap_length = Int(round(window_length / 2)) # audioflux
# # mel
# window_norm = true
# mel_style = :htk
# filterbank_design_domain = :linear
# filterbank_normalization = :bandwidth
# spectrum_type = :power
# num_bands = 32
# frequency_range = [0, sr/2]
# #mfcc
# num_coeffs = 13
# rectification = :log
# log_energy_pos = :append
# delta_window_length = 9

# # set default values for :htk of :slaney, assuming sr = 16000
# if (mel_style == :default_htk)
#     num_bands = 24
#     frequency_range = [0, sr/2]
# elseif (mel_style == :default_slaney)
#     num_bands = 40
#     frequency_range = [133, 6854]
# end

# # setup and data structures definition
# setup = signal_setup(
#     sr=sr,

#     # fft
#     window_type=window_type,
#     window_length=window_length,
#     overlap_length=overlap_length,
#     window_norm=window_norm,

#     # spectrum
#     frequency_range=frequency_range,
#     spectrum_type=spectrum_type,

#     # mel
#     mel_style=mel_style,
#     num_bands=num_bands,
#     filterbank_design_domain=filterbank_design_domain,
#     filterbank_normalization=filterbank_normalization,

#     # mfcc
#     num_coeffs=num_coeffs,
#     rectification=rectification,
#     log_energy_pos=log_energy_pos,
#     delta_window_length=delta_window_length,
# )

# data = signal_data( # no abs, no pre emphasis
#     x=Float64.(x)
# )

# takeFFT(data, setup)
# lin_spectrogram(data, setup)
# mel_spectrogram(data, setup)
# _mfcc(data, setup)
# spectral_features(data, setup)
# f0(data, setup)