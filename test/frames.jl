using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = Audio911.load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn Frames(audiofile)

win = Audio911.movingwindow(winsize=256, winstep=128)
@test_nowarn Frames(audiofile; win, type=hamming, periodic=true)

frames_data = Frames(audiofile; win, type=hamming)
@test length(frames_data) == 132

audiofile = load(wav_file; sr=6000)
frames    = Frames(audiofile)
@test_nowarn get_data(frames)
@test get_size(frames)    == 256
@test get_step(frames)    == 128
@test get_overlap(frames) == 128

audiofile = load(wav_file; sr=10000)
frames    = Frames(audiofile)
@test_nowarn get_data(frames)
@test get_size(frames)    == 512
@test get_step(frames)    == 256
@test get_overlap(frames) == 256

audiofile = Audio911.load(mp3_file; mono=false)
frames    = Frames(audiofile)
@test_nowarn get_data(frames)

# automatic mono conversion
a1 = load(mp3_file; mono=false)
a2 = load(mp3_file; mono=true)
f1 = Frames(a1)
f2 = Frames(a2)
@test get_data(f1) == get_data(f2)

