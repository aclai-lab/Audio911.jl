using Audio911

test_files_dir() = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file     = test_file("test.wav")
mp3_file     = test_file("test.mp3")

audiofile = load(wav_file; mono=true, sr=8000, norm=false)

mel_spec = get_mel_spec(audiofile, mel_style=:htk)
# 2-element Vector{Float64}:
#     0.0
#  2146.06452750619

# 28-element Vector{Float64}:
#     0.0
#    51.1517145741367
#   106.04128329666476
#   164.94184566546735
#   228.14650054076301
#     ⋮
#  2844.69902982251
#  3103.7239341435093
#  3381.676792712254
#  3679.940744547532
#  3999.9999999999995

mel_spec = get_mel_spec(audiofile, mel_style=:slaney)
# 2-element Vector{Float64}:
#   0.0
#  35.163760314616646

# 28-element Vector{Float64}:
#     0.0
#    86.82409954226333
#   173.64819908452665
#   260.47229862678995
#   347.2963981690533
#     ⋮
#  2795.8486475522436
#  3057.7377911682142
#  3344.1582782830383
#  3657.4079773976755
#  4000.0

mel_spec = get_mel_spec(audiofile)