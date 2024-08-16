using Pkg
Pkg.activate("/home/paso/results")
using Revise, Audio911, BenchmarkTools

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

audio1 = audio911features(wavfile)
audio2 = audio911features(wavfile, setup(sr=32000))
analysis1 = audio911features(wavfile, :mfcc, setup(stft_length=1024, mel_bands=52, mfcc_source=:cwt))
analysis2 = audio911features(wavfile, (:cwt, :melfb, :mel), setup(stft_length=1024, mel_bands=52, mfcc_source=:cwt))

# abstracttrees
# unique(collect(AbstractTrees.PreOrderDST(config))
# [i in AbstractTrees.childer(j) for i in nodes, j in nodes)]