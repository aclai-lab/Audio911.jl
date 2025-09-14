using Audio911

@test_nowarn File{format"WAV"}(wav_file)
@test_nowarn File{format"MP3"}(mp3_file)

@test_nowarn load(wav_file)
@test_nowarn load(mp3_file)
# @test_nowarn load(ogg_file)
# @test_nowarn load(flac_file)

audiofile = load(mp3_file)
@test audiofile isa AudioFile
@test nchannels(audiofile) == 1
@test samplerate(audiofile) == 44100

audiofile = load(mp3_file; mono=false)
@test audiofile isa AudioFile
@test nchannels(audiofile) == 2
@test samplerate(audiofile) == 44100

audiofile = load(mp3_file; sr=8000)
@test audiofile isa AudioFile
@test samplerate(audiofile) == 8000

audio      = load(mp3_file; norm=false)
audio_norm = load(mp3_file; norm=true)
@test sum(abs.(data(audio_norm))) > sum(abs.(data(audio)))

@test eltype(audio) == Float32
