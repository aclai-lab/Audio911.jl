using Revise, Audio911, 

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

# ---------------------------------------------------------------------------- #
#      continuous wavelets transform to feed mfcc and spectral features        #
# ---------------------------------------------------------------------------- #
# load audiofile
# sample rate suggested for vocal analysis is 8000hz
# always good pratice to normalize the audio beforehand
audio = load_audio(
    file=wavfile, 
    sr=8000, 
    norm=true,
);

# *** bear in mind that you can always display a plot of waht's going on *** #
# simply by typing "display(plot)" after the function call, like this:
# display(audio)

# if needed you can optionally perform vocal speech detection to cut silence and/or noise
audio = speech_detector(audio=audio);

# compute the continuous wavelet transform filterbank
# as a starting point we suggest to use the bump wavelet
# the frequency range set starting at 80hz is due to wavelets filterbank algorithm:
# starting from 0hz, the filterbank is too compressed at low frequencies.
# our suggest is to avoid the low frequencies, and start at 80hz, for computational efficiency.
cwtfb = get_cwtfb(
    audio=audio,
    wavelet = :bump,
	vpo = 10,
    freq_range = (80, round(Int, audio.sr / 2))
);

# compute the continuous wavelet transform 
# for further mel analysis we suggest to use a power normalization.
# we don't suggest to use a logarithmic scale, as it's intended to be used for visualization.
cwt = get_cwt(
    audio=audio, 
    fbank=cwtfb,
    norm = :power, # :power, :magnitude, :pow2mag
    db_scale = false,
);

# now we are in frequency domain, and we already applied a filterbank, just like the mel spectrogram.
# but bear in mind that we are used to use a number of coefficients that is half of the number of mel bands,
# so we suggest to use the same proportion:
ncoeffs = round(Int, size(cwt.spec, 1) / 2)
# anyway you can skip that because if you don't declare ncoeffs,
# it will be defaulted to half of the number of filterbank coefficients.
mfcc = get_mfcc(
    source=cwt,
    ncoeffs=ncoeffs,
    rectification=:log, # :log, :cubic_root
    dither=true,
);

# now we should continue with fundamental frequency and spectral features
# because if we want a feature matrix, we need to have the same length for every feature.
# TODO resolve current hack in f0 algorithm.
cwtf0 = get_stft(
    audio=audio, 
    stft_length=256,
	win_type = (:hann, :periodic),
	win_length=256,
	overlap_length = 255,
	norm = :power, # :none, :power, :magnitude, :pow2mag
);

f0 = get_f0(
    source=cwtf0,
    method = :nfc,
	freq_range = (80, 400)
);
#[hack]
m = median(f0.f0)
f0.f0 = vcat(f0.f0, zeros(length(audio.data)-length(f0.f0)))

spect = get_spectrals(
    source=cwt,
    freq_range = (0, round(Int, audio.sr / 2))
);

# collect the features in a single Matrix
# case 1: we want rows as features and columns as frames
audio_features_rows = vcat(
    cwt.spec,
    mfcc.mfcc,
    f0.f0',
    spect.centroid',
    spect.crest',
    spect.entropy',
    spect.flatness',
    spect.flux',
    spect.kurtosis',
    spect.rolloff',
    spect.skewness',
    spect.decrease',
    spect.slope',
    spect.spread'
);

# case 2:we want colums as features and rows as frames
audio_features_cols = hcat(
    cwt.spec',
    mfcc.mfcc',
    f0.f0,
    spect.centroid,
    spect.crest,
    spect.entropy,
    spect.flatness,
    spect.flux,
    spect.kurtosis,
    spect.rolloff,
    spect.skewness,
    spect.decrease,
    spect.slope,
    spect.spread
);