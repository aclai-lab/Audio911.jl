using Audio911

test_files_dir() = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file     = test_file("test.wav")
mp3_file     = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

mel_spec = get_melspec(audiofile, mel_style=:htk)
mel_spec = get_melspec(audiofile, mel_style=:slaney)

audiofile = load(wav_file)
frames = get_frames(audiofile, win=MovingWindow(window_size=480, window_step=160), type=(:hamming, :periodic))
stft = get_stft(frames)
mel_spec = get_melspec(stft, mel_bands=32)

mel_spec = get_melspec(audiofile, fb_norm=:area)
mel_spec = get_melspec(audiofile, fb_norm=:none)

