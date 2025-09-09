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

