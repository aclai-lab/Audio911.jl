using Test
using Audio911

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")

# Load audio

sr_resample = 8000
x, sr = load_audio("$TESTPATH/common_voice_en_23616312.wav"; sr = sr_resample)

# No resampling (use original sampling)
# x, sr = load_audio("$TESTPATH/common_voice_en_23616312.wav")



fft_length = 256
frequency_range=Int[0, sr/2]
mel_bands = 26
num_coeffs = 13

setup = FeatureSetup(
    sr=sr,
    # fft
    window_type=[:hann, :periodic],
    window_length=fft_length,
    overlap_length=Int(round(fft_length * 0.500)),
    window_norm=false,
    # spectrum
    frequency_range=frequency_range,
    spectrum_type=:power, # :power, :magnitude
    # mel
    mel_style=:htk, # :htk, :slaney
    mel_bands=mel_bands,
    filterbank_design_domain=:linear,
    filterbank_normalization=:bandwidth, # :bandwidth, :area, :none
    frequency_scale=:mel,
    # mfcc
    num_coeffs=num_coeffs,
    normalization_type=:dithered, # :standard, :dithered
    rectification=:log,
    log_energy_source=:standard, # :standard (after windowing), :mfcc
    log_energy_pos=:none, #:append, :replace, :none
    delta_window_length=9,
    delta_matrix=:standard, # :standard, :transposed
    # spectral
    spectral_spectrum=:linear # :linear, :linear_focused, :mel
)

# convert to Float64
x = Float64.(x)

################################################################################
# Pipelined feature extraction
################################################################################

@test_nowarn data = extractfeatures(x, setup)

features = mel_spectrogram(data)
features = _mfcc(data)
features = lin_spectrogram(data)
features = spectral_features(data)
features = f0(data) # pay attention to fft length

################################################################################
# Standalone feature extraction
################################################################################

features = mel_spectrogram(x)
features = mel_spectrogram(x, setup)
features = mel_spectrogram(x, mel_style = :htk, mel_bands = mel_bands)
features = mel_spectrogram(x, setup, mel_style = :htk, mel_bands = mel_bands)

features = lin_spectrogram(x)

################################################################################

data = extractFFT(setup, x)


data = signal_data(
    x = x,
)

@test_nowarn extractfeatures(data, setup) isa Vector
