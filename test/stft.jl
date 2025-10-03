using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn Stft(audiofile)

@test_nowarn Stft(audiofile; win=MovingWindow(size=1024, step=512))
@test_nowarn Stft(audiofile; win=MovingWindow(size=512, step=256))

win = MovingWindow(size=256, step=128)
frames = AudioFrames(audiofile; win, type=hamming, periodic=true)

@test_nowarn Stft(frames)

@test_nowarn Stft(frames; stft_size=1024)
@test_nowarn Stft(frames; stft_size=2048)
@test_throws ArgumentError Stft(frames; stft_size=128)

@test_nowarn Stft(frames; spectrum_type=power)
@test_nowarn Stft(frames; spectrum_type=magnitude)