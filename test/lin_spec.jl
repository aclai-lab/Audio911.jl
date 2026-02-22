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
matlab_files_dir()    = joinpath(dirname(@__FILE__), "matlab_files/lin_spec")
matlab_file(filename) = joinpath(matlab_files_dir(), filename)

# [audio_wav,fs_wav] = audioread("/home/paso/Documents/Aclai/PasoStudio73/Audio911.jl/test/test_files/test.wav");

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     linearSpectrum=true);

# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audioIn);

matfile = matlab_file("matlab_linspec_01.mat")
mat = MAT.matread(matfile)
mat_lin_spec = mat["features"]

frames = Frames(audiofile; winsize=512, winstep=256, type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
lin_spec = LinSpec(stft; freqrange=(100,1000), win_norm=true)

@test isapprox(get_data(lin_spec), mat_lin_spec)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=bartlett(256), ...
#     OverlapLength=128, ...
#     linearSpectrum=true);

# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[200,800], SpectrumType="magnitude", WindowNormalization=false)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_linspec_02.mat")
mat = MAT.matread(matfile)
mat_lin_spec = mat["features"]

frames = Frames(audiofile; winsize=256, winstep=128, type=bartlett, periodic=false)
stft = Stft(frames; spectrum=magnitude)
lin_spec = LinSpec(stft; freqrange=(200,800), win_norm=false)

@test isapprox(get_data(lin_spec), mat_lin_spec)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hann(417,"periodic"), ...
#     OverlapLength=111, ...
#     linearSpectrum=true);

# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[97,1243], SpectrumType="power", WindowNormalization=false)
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_linspec_03.mat")
mat = MAT.matread(matfile)
mat_lin_spec = mat["features"]

frames = Frames(audiofile; winsize=417, winstep=417-111, type=hanning, periodic=false)
stft = Stft(frames; spectrum=power)
lin_spec = LinSpec(stft; freqrange=(97,1243), win_norm=false)

@test isapprox(get_data(lin_spec), mat_lin_spec)
