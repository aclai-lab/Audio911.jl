using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

winfunc = MovingWindow(size=256, step=128)
frames  = AudioFrames(audiofile; winfunc, type=hamming, periodic=true)