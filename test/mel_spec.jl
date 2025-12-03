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
matlab_files_dir()    = joinpath(dirname(@__FILE__), "matlab_files/mel_spec")
matlab_file(filename) = joinpath(matlab_files_dir(), filename)

# [audio_wav,fs_wav] = audioread("/home/paso/Documents/Aclai/PasoStudio73/Audio911.jl/test/test_files/test.wav");

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     melSpectrum=true);
# setExtractorParameters(aFE,"melSpectrum",FrequencyRange=[100,1000], ... 
#     SpectrumType="power", NumBands=26, FilterBankNormalization="bandwidth", ...
#     WindowNormalization=true, FilterBankDesignDomain="linear", ...
#     MelStyle="oshaughnessy", ApplyLog=false)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_melspec_01.mat")
mat = MAT.matread(matfile)
mat_mel_spec = mat["features"]

frames = AudioFrames(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
mel_spec = MelSpec(stft; win_norm=true, freqrange=(100,1000), nbands=26, norm=bandwidth, domain=:linear, scale=htk)

@test isapprox(get_data(mel_spec), mat_mel_spec)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     melSpectrum=true);
# setExtractorParameters(aFE,"melSpectrum",FrequencyRange=[100,1000], ... 
#     SpectrumType="magnitude", NumBands=26, FilterBankNormalization="bandwidth", ...
#     WindowNormalization=true, FilterBankDesignDomain="linear", ...
#     MelStyle="oshaughnessy", ApplyLog=false)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_melspec_02.mat")
mat = MAT.matread(matfile)
mat_mel_spec = mat["features"]

frames = AudioFrames(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true)
stft = Stft(frames; spectrum=magnitude)
mel_spec = MelSpec(stft; win_norm=true, freqrange=(100,1000), nbands=26, norm=bandwidth, domain=:linear, scale=htk)

@test isapprox(get_data(mel_spec), mat_mel_spec)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     melSpectrum=true);
# setExtractorParameters(aFE,"melSpectrum",FrequencyRange=[100,1000], ... 
#     SpectrumType="power", NumBands=26, FilterBankNormalization="bandwidth", ...
#     WindowNormalization=true, FilterBankDesignDomain="linear", ...
#     MelStyle="oshaughnessy", ApplyLog=false)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_melspec_03.mat")
mat = MAT.matread(matfile)
mat_mel_spec = mat["features"]

frames = AudioFrames(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
mel_spec = MelSpec(stft; win_norm=false, freqrange=(100,1000), nbands=26, norm=bandwidth, domain=:linear, scale=htk)

@test isapprox(get_data(mel_spec), mat_mel_spec)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     melSpectrum=true);
# setExtractorParameters(aFE,"melSpectrum",FrequencyRange=[100,1000], ... 
#     SpectrumType="magnitude", NumBands=26, FilterBankNormalization="bandwidth", ...
#     WindowNormalization=true, FilterBankDesignDomain="linear", ...
#     MelStyle="oshaughnessy", ApplyLog=false)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_melspec_04.mat")
mat = MAT.matread(matfile)
mat_mel_spec = mat["features"]

frames = AudioFrames(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true)
stft = Stft(frames; spectrum=magnitude)
mel_spec = MelSpec(stft; win_norm=false, freqrange=(100,1000), nbands=26, norm=bandwidth, domain=:linear, scale=htk)

@test isapprox(get_data(mel_spec), mat_mel_spec)
