using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = Audio911.load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn Stft(audiofile)

@test_nowarn Stft(audiofile; win=Audio911.movingwindow(winsize=1024, winstep=512))
@test_nowarn Stft(audiofile; win=Audio911.movingwindow(winsize=512, winstep=256))

win = Audio911.movingwindow(winsize=512, winstep=256)
frames = Frames(audiofile; win, type=hamming, periodic=true)

@test_nowarn Stft(frames)

@test_nowarn Stft(frames; nfft=1024)
@test_nowarn Stft(frames; nfft=2048)
@test_throws ArgumentError Stft(frames; nfft=128)

@test_nowarn Stft(frames; spectrum=power)
@test_nowarn Stft(frames; spectrum=magnitude)

stft = Stft(frames; spectrum=power)
