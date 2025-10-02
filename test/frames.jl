using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn AudioFrames(audiofile)

win = MovingWindow(size=256, step=128)
@test_nowarn AudioFrames(audiofile; win, type=hamming, periodic=true)

frames_data = AudioFrames(audiofile; win, type=hamming)
@test length(frames_data) == 132

audiofile = load(wav_file; sr=6000)
frames    = AudioFrames(audiofile)
@test_nowarn get_frames(frames)
@test get_size(frames)    == 256
@test get_step(frames)    == 128
@test get_overlap(frames) == 128

audiofile = load(wav_file; sr=10000)
frames    = AudioFrames(audiofile)
@test_nowarn get_frames(frames)
@test get_size(frames)    == 512
@test get_step(frames)    == 256
@test get_overlap(frames) == 256

audiofile = load(mp3_file; mono=false)
frames    = AudioFrames(audiofile)
@test_nowarn get_frames(frames)

# automatic mono conversion
a1 = load(mp3_file; mono=false)
a2 = load(mp3_file; mono=true)
f1 = AudioFrames(a1)
f2 = AudioFrames(a2)
@test get_frames(f1) == get_frames(f2)

