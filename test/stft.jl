using Audio911

test_files_dir() = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file     = test_file("test.wav")
mp3_file     = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

stft = get_stft(audiofile)
@test stft isa Stft
@test get_stft(stft) isa Matrix{Float64}
@test get_stft_freq(stft) isa Vector{Float64}

info = get_info(stft)
@test info.sr == 8000
@test info.stft_size == 256
@test info.win_size == 256
@test info.overlap == 128
@test info.frequency_range == (0, 4000)
@test info.spectrum_type == :power

stft = get_stft(audiofile, win=(MovingWindow(window_size=1024, window_step=256)), type=(:hann,:periodic))
info = get_info(stft)
@test info.stft_size == 1024
@test info.win_size == 1024
@test info.overlap == 768

@test_nowarn get_stft(audiofile, stft_size=512)

# test invalid frequency range - negative start
@test_throws ArgumentError get_stft(audiofile; frequency_range=(-100, 4000))
# test invalid frequency range - start >= end
@test_throws ArgumentError get_stft(audiofile; frequency_range=(2000, 2000))
# test invalid frequency range - end > sr/2
@test_throws ArgumentError get_stft(audiofile; frequency_range=(0, 5000))  # sr/2 = 4000
# test invalid overlap (create frames with zero overlap)
no_overlap_frames = get_frames(audiofile; win=MovingWindow(window_size=256, window_step=256))
@test_throws ArgumentError get_stft(no_overlap_frames, 8000)
# test invalid stft_size < window_size
@test_throws ArgumentError get_stft(audiofile; stft_size=128)  # window_size is 256
# test stereo audio (should fail for STFT)
stereofile = load(mp3_file; mono=false)
@test_throws ArgumentError get_stft(stereofile)
# test invalid spectrum_type
@test_throws AssertionError get_stft(audiofile; spectrum_type=:invalid)
