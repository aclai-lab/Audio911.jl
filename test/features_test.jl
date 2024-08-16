using Pkg
Pkg.activate("/home/paso/results")
using Revise, Audio911, BenchmarkTools

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

# --- audio ------------------------------------------------------------------ #
sr = 16000
audio = load_audio(file=wavfile);
audio = load_audio(file=wavfile, sr=sr, norm=true);
show(audio)
plot(audio)
save_audio(audio=audio, file="/home/paso/Documents/testnorm.wav")

# --- stft ------------------------------------------------------------------- #
stftspec = get_stft(audio=audio);
stftspec = get_stft(audio=audio, stft_norm=:magnitude);
show(stftspec)
plot(stftspec)

linspec = get_linspec(source=stftspec);
linspec = get_linspec(source=stftspec, db_scale=true);
display(linspec)

melfb = get_melfb(stft=stftspec);
melfb = get_melfb(stft=stftspec, scale=:erb);
# display(melfb);

melspec =  get_melspec(stft=stftspec, fbank=melfb);
melspec =  get_melspec(stft=stftspec, fbank=melfb, db_scale=true);
display(melspec)

melspec =  get_melspec(stft=stftspec, fbank=melfb);
mfcc = get_mfcc(source=melspec);
mfcc = get_mfcc(source=melspec, rectification=:cubic_root);
display(mfcc)

deltas = get_deltas(source=mfcc);
deltas = get_deltas(source=mfcc, d_length=7);
display(deltas)

spect = get_spectrals(source=stftspec);
spect = get_spectrals(source=stftspec, freq_range=(100, 2000));
display(spect)

f0 = get_f0(source=stftspec);
f0 = get_f0(source=stftspec, freq_range=(300, 2000));
display(f0)

cwtfb = get_cwtfb(audio=audio);
cwtfb = get_cwtfb(audio=audio, wavelet=:bump, freq_range=(100,4000));
display(cwtfb)

cwt = get_cwt(audio=audio, fbank=cwtfb);
cwt = get_cwt(audio=audio, fbank=cwtfb, norm=:magnitude, db_scale=true);
display(cwt)

# BROKEN!
# cwtspec = get_linspec(source=cwt);
# cwtspec = get_linspec(source=cwt, db_scale=true);
# display(cwtspec)