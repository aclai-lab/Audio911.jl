using Audio911

test_files_dir() = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file  = test_file("test.wav")
mp3_file  = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn get_frames(audiofile)

win = MovingWindow(window_size=256, window_step=128)
type = (:hann, :periodic)

@test_nowarn get_frames(audiofile; win, type)

frames_data = get_frames(audiofile; win, type)
@test length(frames_data) == 132

audiofile = load(wav_file; sr=6000)
frames = get_frames(audiofile)
@test_nowarn get_frames(frames)
@test get_wsize(frames)  == 256
@test get_wstep(frames)  == 128
@test get_ovrlap(frames) == 128

audiofile = load(wav_file; sr=10000)
frames = get_frames(audiofile)
@test_nowarn get_frames(frames)
@test get_wsize(frames)  == 512
@test get_wstep(frames)  == 256
@test get_ovrlap(frames) == 256

audiofile = load(mp3_file; mono=false)
frames = get_frames(audiofile)
@test_nowarn get_frames(frames)
@test nchannels(frames) == 2
