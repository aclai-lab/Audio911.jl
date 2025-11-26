using Test
using Audio911

using MAT

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = Audio911.load(wav_file, format=Float64)

frames = AudioFrames(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true)
stft = Stft(frames; spectrum_type=power)

fbank = FBank(stft; nbands=26, norm=:bandwidth, domain=:linear, freqrange=(100,1000))

@show Threads.nthreads()

@btime FBank(stft)
# 23.573 μs (79 allocations: 60.27 KiB)
