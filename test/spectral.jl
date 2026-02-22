using Test
using Audio911

using MAT

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = Audio911.load(wav_file, format=Float64)
@test get_sr(audiofile) == 16000
@test is_norm(audiofile) == false
@test size(get_data(audiofile), 2) == 1

# ---------------------------------------------------------------------------- #
#                            test against matlab                               #
# ---------------------------------------------------------------------------- #
matlab_files_dir()    = joinpath(dirname(@__FILE__), "matlab_files/spectral")
matlab_file(filename) = joinpath(matlab_files_dir(), filename)

# [audio_wav,fs_wav] = audioread("/home/paso/Documents/Aclai/Audio911.jl/test/test_files/test.wav");

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralCentroid=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralCentroid.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

stft = Stft(audiofile; winsize=512, winstep=256, type=hamming, periodic=true, spectrum=power)
lin_spec = LinSpec(stft; freqrange=(100,1000), win_norm=true)
spectral = SpectralCentroid(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralCrest=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralCrest.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralCrest(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralDecrease=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralDecrease.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralDecrease(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralEntropy=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralEntropy.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralEntropy(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralFlatness=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralFlatness.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralFlatness(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralFlux=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# setExtractorParameters(aFE,"spectralFlux",  NormType=1)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralFlux1.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralFlux(lin_spec; p=1)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralFlux=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# setExtractorParameters(aFE,"spectralFlux",  NormType=2)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralFlux2.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralFlux(lin_spec, p=2)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralKurtosis=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralKurtosis.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralKurtosis(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralRolloffPoint=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# setExtractorParameters(aFE,"spectralRolloffPoint", Threshold=0.95)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralRolloff95.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralRolloff(lin_spec; threshold=0.95)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralRolloffPoint=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# setExtractorParameters(aFE,"spectralRolloffPoint", Threshold=0.65)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralRolloff65.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralRolloff(lin_spec, threshold=0.65)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralSkewness=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralSkewness.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralSkewness(lin_spec)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralSlope=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralSlope.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralSlope(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     spectralSpread=true);
# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_spectralSpread.mat")
mat = MAT.matread(matfile)
mat_spectral = mat["features"]

spectral = SpectralSpread(lin_spec)

@test isapprox(get_data(spectral), mat_spectral)
