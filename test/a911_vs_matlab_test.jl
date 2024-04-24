using JLD2, DataFrames
# using Audio911

#--------------------------------------------------------------------------------------#
#                                 experiment settings                                  #
#--------------------------------------------------------------------------------------#
# ds_type = :gender
# # ds_type = :age2bins

# profile = :full
# # profile = :gender
# # profile = :age

# n_samples = 1000

# # audio parameters
# sr = 8000
# fft_length = 256
# frequency_range = (0, floor(Int, sr/2))
# mel_bands = 26
# num_coeffs = 13


#--------------------------------------------------------------------------------------#
#                                       audio911                                       #
#--------------------------------------------------------------------------------------#
# setup = AudioSetup(
#     sr=sr,
#     # fft
#     window_type=(:hann, :periodic),
#     window_length=fft_length,
#     overlap_length=Int(round(fft_length * 0.500)),
#     window_norm=false,
#     # spectrum
#     frequency_range=frequency_range,
#     spectrum_type=:power, # :power, :magnitude
#     # mel
#     mel_style=:htk, # :htk, :slaney
#     mel_bands=mel_bands,
#     filterbank_design_domain=:linear,
#     filterbank_normalization=:bandwidth, # :bandwidth, :area, :none
#     frequency_scale=:mel,
#     # mfcc
#     num_coeffs=num_coeffs,
#     normalization_type=:standard, # :standard, :dithered
#     rectification=:log,
#     log_energy_source=:standard, # :standard (after windowing), :mfcc
#     log_energy_pos=:none, #:append, :replace, :none
#     delta_window_length=9,
#     delta_matrix=:standard, # :standard, :transposed
#     # spectral
#     spectral_spectrum=:linear # :linear, :mel
# )

# data = AudioData(
#     x=x
# )

# get_fft!(data, setup)
# mel_spectrogram(data, setup)
# _mfcc(data, setup)
# lin_spectrogram(data, setup)
# spectral_features(data, setup)
# f0(data, setup)

# return vcat((
#     data.mel_spectrogram',
#     data.mfcc_coeffs',
#     data.mfcc_delta',
#     data.mfcc_deltadelta',
#     data.spectral_centroid',
#     data.spectral_crest',
#     data.spectral_decrease',
#     data.spectral_entropy',
#     data.spectral_flatness',
#     data.spectral_flux',
#     data.spectral_kurtosis',
#     data.spectral_rolloff',
#     data.spectral_skewness',
#     data.spectral_slope',
#     data.spectral_spread',
#     data.f0'
# )...)

#--------------------------------------------------------------------------------------#
#                                      load jld2                                       #
#--------------------------------------------------------------------------------------#
a911_file = "/home/riccardopasini/Documents/Aclai/Datasets/Common_voice_ds/debug/test 4 spectral/spcds_gender_audio911_mag.jld2"
matlab_file = "/home/riccardopasini/Documents/Aclai/Datasets/Common_voice_ds/debug/test 4 spectral/spcds_gender_audio911_matlab.jld2"

a911_jld2 = jldopen(a911_file)
matlab_jld2 = jldopen(matlab_file)

dfa, ya = a911_jld2["dataframe_validated"]
dfm, ym = matlab_jld2["dataframe_validated"]

#--------------------------------------------------------------------------------------#
#                                        test                                          #
#--------------------------------------------------------------------------------------#
if ya != ym
    error("Y labels are different!")
end

if size(dfa, 1) != size(dfm, 1) || size(dfa, 2) != size(dfm, 2)
    error("Sizes are different!")
end

for i in axes(dfa, 2)
    if !isapprox(dfa[!, i], dfm[!, i])
        @warn("col $i is different!")
        # for j in axes(dfa, 1)
        #     if !isapprox(dfa[j, i], dfm[j, i])
        #         @warn("row $j is different!") 
        #         for k in 1:length(dfa[j,i])
        #             # println(dfa[j, i][k], " vs ", dfm[j, i][k])
        #             if !isapprox(dfa[j, i][k], dfm[j, i][k])
        #                 @warn("data $k is different: ", dfa[j, i][k], " vs ", dfm[j, i][k]) 
        #             end
        #         end
        #     end
        # end
    end
end
