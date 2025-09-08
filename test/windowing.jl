using Audio911

test_files_dir() = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file     = test_file("test.wav")
mp3_file     = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

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