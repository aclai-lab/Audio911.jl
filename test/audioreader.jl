using Audio911

@test_nowarn File{format"WAV"}(wav_file)
@test_nowarn File{format"MP3"}(mp3_file)

@test_nowarn load(wav_file)
@test_nowarn load(mp3_file)
# @test_nowarn load(ogg_file)
# @test_nowarn load(flac_file)