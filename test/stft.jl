using Test
using Audio911

using MAT

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = Audio911.load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn Stft(audiofile)

@test_nowarn Stft(audiofile; win=MovingWindow(size=1024, step=512))
@test_nowarn Stft(audiofile; win=MovingWindow(size=512, step=256))

win = MovingWindow(size=512, step=256)
frames = AudioFrames(audiofile; win, type=hamming, periodic=true)

@test_nowarn Stft(frames)

@test_nowarn Stft(frames; stft_size=1024)
@test_nowarn Stft(frames; stft_size=2048)
@test_throws ArgumentError Stft(frames; stft_size=128)

@test_nowarn Stft(frames; spectrum_type=power)
@test_nowarn Stft(frames; spectrum_type=magnitude)

stft = Stft(frames; spectrum_type=power)

# ---------------------------------------------------------------------------- #
#                            test against matlab                               #
# ---------------------------------------------------------------------------- #
# [audio_wav,fs_wav] = audioread("/home/paso/Documents/Aclai/PasoStudio73/Audio911.jl/test/test_files/test.wav");

# aFE = audioFeatureExtractor( ...
#     SampleRate=fs_wav, ...
#     Window=hamming(512,"periodic"), ...
#     OverlapLength=256, ...
#     linearSpectrum=true);

# setExtractorParameters(aFE,"linearSpectrum",FrequencyRange=[100,1000], SpectrumType="power", WindowNormalization=true)
# features = extract(aFE, audioIn);


audiofile = Audio911.load(wav_file, format=Float64)
@test get_sr(audiofile) == 16000
@test is_norm(audiofile) == false
@test size(get_data(audiofile), 2) == 1

win = MovingWindow(size=512, step=256)
frames = AudioFrames(audiofile; win, type=hamming, periodic=true)