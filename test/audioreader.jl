using Test
using Audio911

using MAT
using Plots

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

# ---------------------------------------------------------------------------- #
#                                audio reader                                  #
# ---------------------------------------------------------------------------- #
@testset "audioreader" begin
    @test_nowarn Audio911.File{format"WAV"}(wav_file)
    @test_nowarn Audio911.File{format"MP3"}(mp3_file)

    @test_nowarn Audio911.load(wav_file)
    @test_nowarn Audio911.load(mp3_file)

    audiofile = Audio911.load(wav_file; format=Float64)
    @test Audio911.get_data(audiofile) isa Array{Float64}
    @test Audio911.get_sr(audiofile) == 16000
    @test Audio911.get_origin_sr(audiofile) == 16000
    @test Audio911.get_nchannels(audiofile) == 1
    @test Audio911.is_norm(audiofile) == false

    audiofile = Audio911.load(mp3_file; norm=true)
    @test Audio911.get_data(audiofile) isa Array{Float32}
    @test Audio911.get_sr(audiofile) == 44100
    @test Audio911.get_origin_sr(audiofile) == 44100
    @test Audio911.get_nchannels(audiofile) == 1
    @test Audio911.is_norm(audiofile) == true

    audiofile = Audio911.load(mp3_file; mono=false)
    @test Audio911.get_data(audiofile) isa Array{Float32}
    @test Audio911.get_sr(audiofile) == 44100
    @test Audio911.get_origin_sr(audiofile) == 44100
    @test Audio911.get_nchannels(audiofile) == 2
    @test Audio911.is_norm(audiofile) == false

    audiofile = Audio911.load(wav_file; sr=48000)
    @test Audio911.get_data(audiofile) isa Array{Float32}
    @test Audio911.get_sr(audiofile) == 48000
    @test Audio911.get_origin_sr(audiofile) == 16000
    @test Audio911.get_nchannels(audiofile) == 1
    @test Audio911.is_norm(audiofile) == false

    audiofile = Audio911.load(mp3_file; mono=false, sr=8000)  
    @test Audio911.get_nchannels(audiofile) == 2 
    @test eltype(audiofile) == Float32
    @test length(audiofile) == 44513
end

# ---------------------------------------------------------------------------- #
#                            test against matlab                               #
# ---------------------------------------------------------------------------- #
matlab_files_dir()    = joinpath(dirname(@__FILE__), "matlab_files/audioread")
matlab_file(filename) = joinpath(matlab_files_dir(), filename)

@testset "against matlab" begin
    matfile_wav = matlab_file("matlab_audioread_wav.mat")
    matfile_mp3 = matlab_file("matlab_audioread_mp3.mat")
    resampled_wav = matlab_file("matlab_resampled_wav.mat")

    mat_wav = MAT.matread(matfile_wav)
    wav_data_mat = mat_wav["audio_wav"]

    a911_wav = Audio911.load(wav_file, format=Float64)
    wav_data_a911 = Audio911.get_data(a911_wav)

    @test isapprox(wav_data_a911, wav_data_mat)

    mat_mp3 = MAT.matread(matfile_mp3)
    mp3_data_mat = mat_mp3["audio_mp3"]

    a911_mp3 = Audio911.load(mp3_file, mono=false)
    mp3_data_a911 = Audio911.get_data(a911_mp3)

    @test isapprox(mp3_data_mat, mp3_data_a911)
end

# ---------------------------------------------------------------------------- #
#                                  resampling                                  #
# ---------------------------------------------------------------------------- #
orig_file = Audio911.load(wav_file)
orig_data = Audio911.get_data(orig_file)
res_file  = Audio911.load(wav_file, sr=8000)
res_data  = Audio911.get_data(res_file)

# plot comparison
plot(orig_data[1:2000], label="Original (16kHz)", title="Comparison", linewidth=1.5)
plot!(1:2:2000, res_data[1:1000], label="Resampled (8kHz)", linewidth=1.5, linestyle=:dash)

orig_file = Audio911.load(mp3_file)
orig_data = Audio911.get_data(orig_file)
res_file  = Audio911.load(mp3_file, sr=96000)
res_data  = Audio911.get_data(res_file)

# calculate the ratio: 96000/44100 ≈ 2.177
ratio = 96000 / 44100

# plot comparison
n_samples = 4000
orig_indices = 1:n_samples
res_indices = round.(Int, 1:ratio:n_samples*ratio)

plot(orig_indices, orig_data[1:n_samples], 
    label="Original (44.1kHz)", 
    title="MP3 Comparison: 44.1kHz vs 96kHz", 
    linewidth=1.5,
    xlabel="Sample index (normalized)",
    ylabel="Amplitude")
plot!(orig_indices, res_data[res_indices], 
    label="Resampled (96kHz)", 
    linewidth=1.5, 
    linestyle=:dash,
    alpha=0.7)
