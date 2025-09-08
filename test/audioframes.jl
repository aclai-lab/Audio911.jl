using Audio911

test_files_dir() = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file     = test_file("test.wav")
mp3_file     = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

@test_nowarn Audio911.get_frames(audiofile)

windows = (
    (AdaptiveWindow(nwindows=3, relative_overlap=0.3), 3),
    (WholeWindow(), 1),
    (SplitWindow(nwindows=2), 2),
    (MovingWindow(window_size=256, window_step=128), 132),
)
type = (:hann, :periodic)

for win in windows
    @test_nowarn Audio911.get_frames(audiofile; win=win[1], type)
end

for win in windows
    frames_data = Audio911.get_frames(audiofile; win=win[1], type)
    @test length(frames_data) == win[2]
end

audiofile = load(wav_file; sr=6000)
frames = Audio911.get_frames(audiofile)
@test frames.win.params.window_size == 256
@test frames.win.params.window_step == 128

audiofile = load(wav_file; sr=10000)
frames = Audio911.get_frames(audiofile)
@test frames.win.params.window_size == 512
@test frames.win.params.window_step == 256