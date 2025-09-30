using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = load(wav_file)
winfunc   = MovingWindow(size=512, step=256)

audioframes = get_frames(
	audiofile;
    winfunc,
	type=hamming,
    periodic=true
)
