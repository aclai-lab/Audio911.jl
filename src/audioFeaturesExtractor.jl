function audio_features_extractor(
    x::AbstractVector{T};
    sr::Int64,
    dataset::Symbol=:all,

    # fft
    fft_length::Int64=256,
    window_type::Vector{Symbol}=[:hann, :periodic],
    window_length::Int64=fft_length,
    overlap_length::Int64=Int(round(fft_length * 0.500)),
    # window_length::Int64 = Int(round(0.03 * sr)),
    # overlap_length::Int64 = Int(round(0.02 * sr)),
    window_norm::Bool=:false,

    # spectrum
    frequency_range::Vector{Int64}=Int[0, sr/2],
    spectrum_type::Symbol=:power,

    # mel
    mel_style::Symbol=:htk,
    mel_bands::Int64=26,
    filterbank_design_domain::Symbol=:linear,
    filterbank_normalization::Symbol=:bandwidth,
    frequency_scale::Symbol=:mel,

    # mfcc
    num_coeffs::Int64=13,
    normalization_type::Symbol=:dithered,
    rectification::Symbol=:log,
    log_energy_source::Symbol=:standard,
    log_energy_pos::Symbol=:replace,
    delta_window_length::Int64=9,
    delta_matrix::Symbol=:transposed,

    # spectral
    spectral_spectrum::Symbol=:linear
) where {T<:AbstractFloat}

    setup = signal_setup(
        sr=sr,

        # fft
        window_type=window_type,
        window_length=window_length,
        overlap_length=overlap_length,
        window_norm=window_norm,

        # spectrum
        frequency_range=frequency_range,
        spectrum_type=spectrum_type,

        # mel
        mel_style=mel_style,
        mel_bands=mel_bands,
        filterbank_design_domain=filterbank_design_domain,
        filterbank_normalization=filterbank_normalization,
        frequency_scale=frequency_scale,

        # mfcc
        num_coeffs=num_coeffs,
        normalization_type=normalization_type,
        rectification=rectification,
        log_energy_source=log_energy_source,
        log_energy_pos=log_energy_pos,
        delta_window_length=delta_window_length,
        delta_matrix=delta_matrix,

        # spectral
        spectral_spectrum=spectral_spectrum
    )

    # convert to Float64
    x = Float64.(x)

    # preemphasis
    # not siutable for our kind of experiments, maybe for speaker recognition: needs to look over it.
    # zi = 2 * x[1] - x[2]
    # filt!(x, [1.0, -0.97], 1.0, x, [zi])
    # normalize
    # x = x ./ maximum(abs.(x))

    data = signal_data(
        x=x
    )

    takeFFT(data, setup)
    lin_spectrogram(data, setup)
    mel_spectrogram(data, setup)
    _mfcc(data, setup)
    spectral_features(data, setup)
    f0(data, setup) # pay attention to fft length!

    # TODO verificare che il sample sia di lunghezza superiore a fft_length

    if dataset == :full
        vcat((
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
            data.spectral_spread', 
            data.f0'
        )...)

    elseif dataset == :gender
        vcat((
            data.mfcc_coeffs',
            data.mfcc_delta',
            data.mel_spectrogram[:, 1:13]'
        )...)

    elseif dataset == :speaker_recognition
        vcat((
            # data.mel_spectrogram[:, 1:26]',
            data.mfcc_coeffs',
            data.mfcc_delta',
            data.mfcc_deltadelta',
            # data.spectral_centroid', #NO
            data.spectral_crest', #FORSE
            data.spectral_decrease', #SI
            # data.spectral_entropy', #NO
            # data.spectral_flatness', #NO
            data.spectral_flux', #FORSE
            data.spectral_kurtosis', #FORSE
            # data.spectral_rolloff', #NO
            data.spectral_skewness', #FORSE
            # data.spectral_slope', #NO
            # data.spectral_spread', #NO
            data.f0' #SI
        )...)

    else
        error("Unknown dataset type: $dataset.")
    end

end

# bitmask approach
# 1 - mel spectrogram    # every audio parameter is defaulted to optimized value
# 3 - linear spectrogram
# 4 - mfcc
# 5 - delta
# 6 - delta delta
# 7 - centroid
# 8 - crest
# 9 - ecrease
# 10 - entropy
# 11 - flatness
# 12 - flux
# 13 - kurtosis
# 14 - rolloff
# 15 - skewness
# 16 - slope
# 17 - spread
# 18 - f0

# :full = features_bitmask = UInt8[
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1,
#     1
# ]

# function audio_features_extractor(
#     x::AbstractVector{T};
#     sr::Int64,
#     features_bitmask::Vector{UInt8}=UInt8[
#         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1,
#         1
#     ],

#     # fft
#     fft_length::Int64 = 256,
#     window_type::Vector{Symbol}=[:hann, :periodic],
#     window_length::Int64=fft_length,
#     overlap_length::Int64=Int(round(fft_length * 0.500)),
#     # window_length::Int64 = Int(round(0.03 * sr)),
#     # overlap_length::Int64 = Int(round(0.02 * sr)),
#     window_norm::Bool=:false,

#     # spectrum
#     frequency_range::Vector{Int64}=Int[0, sr/2],
#     spectrum_type::Symbol=:power,

#     # mel
#     mel_style::Symbol=:htk,
#     mel_bands::Int64=26,
#     filterbank_design_domain::Symbol=:linear,
#     filterbank_normalization::Symbol=:bandwidth,
#     frequency_scale::Symbol=:mel,

#     # mfcc
#     num_coeffs::Int64=13,
#     normalization_type::Symbol=:dithered,
#     rectification::Symbol=:log,
#     log_energy_source::Symbol=:standard,
#     log_energy_pos::Symbol=:replace,
#     delta_window_length::Int64=9,
#     delta_matrix::Symbol=:transposed,

#     # spectral
#     spectral_spectrum::Symbol=:linear
# ) where {T<:AbstractFloat}

#     # function get_fft_length(sr::Int64)
#     #     sr <= 4000 && return 128
#     #     sr <= 8000 && return 256
#     #     sr <= 16000 && return 512
#     #     return 1024
#     # end

#     # fft_length = get_fft_length(sr)

#     setup = signal_setup(
#         sr=sr,

#         # fft
#         window_type=window_type,
#         window_length=window_length,
#         overlap_length=overlap_length,
#         window_norm=window_norm,

#         # spectrum
#         frequency_range=frequency_range,
#         spectrum_type=spectrum_type,

#         # mel
#         mel_style=mel_style,
#         mel_bands=mel_bands,
#         filterbank_design_domain=filterbank_design_domain,
#         filterbank_normalization=filterbank_normalization,
#         frequency_scale=frequency_scale,

#         # mfcc
#         num_coeffs=num_coeffs,
#         normalization_type=normalization_type,
#         rectification=rectification,
#         log_energy_source=log_energy_source,
#         log_energy_pos=log_energy_pos,
#         delta_window_length=delta_window_length,
#         delta_matrix=delta_matrix,

#         # spectral
#         spectral_spectrum=spectral_spectrum
#     )

#     # convert to Float64
#     x = Float64.(x)

#     # preemphasis
#     # not siutable for our kind of experiments, maybe for speaker recognition: needs to look over it.
#     # zi = 2 * x[1] - x[2]
#     # filt!(x, [1.0, -0.97], 1.0, x, [zi])
#     # normalize
#     # x = x ./ maximum(abs.(x))

#     data = signal_data(
#         x=x
#     )

#     takeFFT(data, setup)
#     lin_spectrogram(data, setup)
#     mel_spectrogram(data, setup)
#     _mfcc(data, setup)
#     spectral_features(data, setup)
#     f0(data, setup) # pay attention to fft length!

#     full_feats = vcat((
#         data.mel_spectrogram',
#         data.mfcc_coeffs',
#         data.mfcc_delta',
#         data.mfcc_deltadelta',
#         data.spectral_centroid',
#         data.spectral_crest',
#         data.spectral_decrease',
#         data.spectral_entropy',
#         data.spectral_flatness',
#         data.spectral_flux',
#         data.spectral_kurtosis',
#         data.spectral_rolloff',
#         data.spectral_skewness',
#         data.spectral_slope',
#         data.spectral_spread', data.f0'
#     )...)
# end