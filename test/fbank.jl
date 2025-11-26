using Test
using Audio911

using MAT

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

matlab_files_dir()    = joinpath(dirname(@__FILE__), "matlab_files/fbanks")
matlab_file(filename) = joinpath(matlab_files_dir(), filename)

# [filterBank, Fc, BW] = designAuditoryFilterBank(fs_wav, ...
#     FrequencyScale="mel", FFTLength=1024, ...
#     NumBands=26, FrequencyRange=[100,1000], ...
#     Normalization="bandwidth", FilterBankDesignDomain="linear", ...
#     MelStyle="oshaughnessy")
# save fb01.mat filterBank Fc BW

fbank = FBank(16000; nfft=1024, nbands=26, norm=bandwidth, domain=:linear, freqrange=(100,1000))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

matfile = matlab_file("fb01.mat")
mat = MAT.matread(matfile)
fbm, fm, bwm = mat["filterBank"], mat["Fc"], mat["BW"]

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

audiofile = Audio911.load(wav_file, format=Float64)
frames = AudioFrames(audiofile; win=movingwindow(winsize=1024, winstep=512), type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
fbank = FBank(stft; nbands=26, norm=bandwidth, domain=:linear, freqrange=(100,1000))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

# [filterBank, Fc, BW] = designAuditoryFilterBank(fs_wav, ...
#     FrequencyScale="mel", FFTLength=733, ...
#     NumBands=15, FrequencyRange=[53,1743], ...
#     Normalization="area", FilterBankDesignDomain="linear", ...
#     MelStyle="slaney")
# save fb02.mat filterBank Fc BW

fbank = FBank(16000; nfft=733, nbands=15, scale=:slaney, norm=area, domain=:linear, freqrange=(53,1743))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

matfile = matlab_file("fb02.mat")
mat = MAT.matread(matfile)
fbm, fm, bwm = mat["filterBank"], mat["Fc"], mat["BW"]

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

audiofile = Audio911.load(wav_file, format=Float64)
frames = AudioFrames(audiofile; win=movingwindow(winsize=733, winstep=733÷2), type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
fbank = FBank(stft; nbands=15, scale=:slaney, norm=area, domain=:linear, freqrange=(53,1743))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

# [filterBank, Fc, BW] = designAuditoryFilterBank(fs_wav, ...
#     FrequencyScale="mel", FFTLength=1024, ...
#     NumBands=26, FrequencyRange=[100,1000], ...
#     Normalization="bandwidth", FilterBankDesignDomain="warped", ...
#     MelStyle="oshaughnessy")
# save fb03.mat filterBank Fc BW

fbank = FBank(16000; nfft=1024, nbands=26, scale=:htk, norm=bandwidth, domain=:warped, freqrange=(100,1000))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

matfile = matlab_file("fb03.mat")
mat = MAT.matread(matfile)
fbm, fm, bwm = mat["filterBank"], mat["Fc"], mat["BW"]

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

audiofile = Audio911.load(wav_file, format=Float64)
frames = AudioFrames(audiofile; win=movingwindow(winsize=1024, winstep=512), type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
fbank = FBank(stft; nbands=26, norm=bandwidth, domain=:warped, freqrange=(100,1000))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

# [filterBank, Fc, BW] = designAuditoryFilterBank(fs_wav, ...
#     FrequencyScale="erb", FFTLength=733, ...
#     NumBands=15, FrequencyRange=[53,1743], ...
#     Normalization="area", FilterBankDesignDomain="linear", ...
#     MelStyle="slaney")
# save fb04.mat filterBank Fc BW

fbank = FBank(16000; nfft=733, nbands=15, scale=:erb, norm=area, freqrange=(53,1743))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

matfile = matlab_file("fb04.mat")
mat = MAT.matread(matfile)
fbm, fm, bwm = mat["filterBank"], mat["Fc"], mat["BW"]

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))

audiofile = Audio911.load(wav_file, format=Float64)
frames = AudioFrames(audiofile; win=movingwindow(winsize=733, winstep=733÷2), type=hamming, periodic=true)
stft = Stft(frames; spectrum=power)
fbank = FBank(stft; nbands=15, scale=:erb, norm=area, freqrange=(53,1743))
fb, f, bw = get_data(fbank), get_freq(fbank), get_bandwidth(fbank)

@test isapprox(fb, fbm)
@test isapprox(f, vec(fm))
@test isapprox(bw, vec(bwm))