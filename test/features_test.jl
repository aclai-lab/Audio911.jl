using Revise, Audio911, BenchmarkTools

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
# abstracttrees
# unique(collect(AbstractTrees.PreOrderDST(config))
# [i in AbstractTrees.childer(j) for i in nodes, j in nodes)]
wavfile = joinpath(TESTPATH, TESTFILE)

sr = 16000
audio = load_audio(file=wavfile);
audio = load_audio(file=wavfile, sr=sr);
display(audio)

stftspec = get_stft(audio=audio);
stftspec = get_stft(audio=audio, norm=:magnitude);
display(stftspec)

linspec = get_linspec(source=stftspec);
linspec = get_linspec(source=stftspec, db_scale=true);
display(linspec)

melfb = get_melfb(stft=stftspec);
melfb = get_melfb(stft=stftspec, scale=:erb);
display(melfb)

melspec =  get_melspec(stft=stftspec, fbank=melfb);
melspec =  get_melspec(stft=stftspec, fbank=melfb, db_scale=true);
display(melspec)

melspec =  get_melspec(stft=stftspec, fbank=melfb);
mfcc = get_mfcc(melspec=melspec);
mfcc = get_mfcc(melspec=melspec, rectification=:cubic_root);
display(mfcc)

deltas = get_deltas(mfcc=mfcc);
deltas = get_deltas(mfcc=mfcc, d_length=7);
display(deltas)

spect = get_spectrals(stft=stftspec);
spect = get_spectrals(stft=stftspec, freq_range=(100, 2000));
display(spect)

f0 = get_f0(stft=stftspec);
f0 = get_f0(stft=stftspec, freq_range=(300, 2000));
display(f0)

cwtfb = get_cwtfb(audio=audio);
cwtfb = get_cwtfb(audio=audio, wavelet=:bump, freq_range=(100,4000));
display(cwtfb)

cwt = get_cwt(audio=audio, fbank=cwtfb);
cwt = get_cwt(audio=audio, fbank=cwtfb, norm=:magnitude, db_scale=true);
display(cwt)

cwtspec = get_linspec(source=cwt);
cwtspec = get_linspec(source=cwt, db_scale=true);
display(cwtspec)