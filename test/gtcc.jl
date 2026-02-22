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
matlab_files_dir()    = joinpath(dirname(@__FILE__), "matlab_files/gtcc")
matlab_file(filename) = joinpath(matlab_files_dir(), filename)

# [audio_wav,fs_wav] = audioread("/home/paso/Documents/Aclai/PasoStudio73/Audio911.jl/test/test_files/test.wav");

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     gtcc=true);
# setExtractorParameters(aFE,"gtcc", NumCoeffs=13, Rectification="log");
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_gtcc_01.mat")
mat = MAT.matread(matfile)
mat_gtcc = mat["features"]

frames = Frames(audiofile; winsize=512, winstep=256, type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
erb_spec = ErbSpec(stft; win_norm=true, nbands=34, norm=bandwidth)
gtcc = Gtcc(erb_spec; ncoeffs=13, rect=mlog)

@test isapprox(get_data(gtcc), mat_gtcc)

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     gtcc=true);
# setExtractorParameters(aFE,"gtcc", NumCoeffs=30, Rectification="cubic-root");
# features = extract(aFE, audio_wav);

matfile = matlab_file("matlab_gtcc_02.mat")
mat = MAT.matread(matfile)
mat_gtcc = mat["features"]

gtcc = Gtcc(erb_spec; ncoeffs=30, rect=cubic_root)

@test isapprox(get_data(gtcc), mat_gtcc)
