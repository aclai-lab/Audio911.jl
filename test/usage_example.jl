using Revise, Audio911

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

# ---------------------------------------------------------------------------- #
#                                 itadata 2024                                 #
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

# compute stft spectrogram
# for vocal analysis we suggest to use a stft length of 256 samples @8000hz and 512 samples @16000hz.
# we also suggest to use a hann, periodic window,
# and using a window of the same length as the stft length, to avoid edge effects.
# overlap should be half of the window length, according to Hanning recommendation.
# for further mel analysis we suggest to use a power normalization.
stft_length	= audio.sr <= 8000 ? 256 : 512

stftspec = get_stft(
    audio=audio, 
    stft_length=stft_length,
	win_type = (:hann, :periodic),
	win_length=stft_length,
	overlap_length = round(Int, stft_length / 2),
	norm = :power, # :power, :magnitude, :pow2mag
);

# compute mel spectrogram part 1: create the filterbank
# we use a filterbank composed of 26 bands, using the mel scale and the htk formula.
# we also suggest to use a bandwidth normalization, and broadband frequency range.
melfb = get_melfb(
    stft=stftspec,
    nbands = 26,
    scale = :mel_htk, # :mel_htk, :mel_slaney, :erb, :bark
    norm = :bandwidth, # :bandwidth, :area, :none
    freq_range = (0, round(Int, audio.sr / 2))
);

# part 2: compute the mel spectrogram
# we don't suggest to use a logarithmic scale, as it's intended to be used for visualization.
melspec =  get_melspec(
    stft=stftspec, 
    fbank=melfb,
    db_scale = false
);

# compute mfcc
# with a mel spectrogram composed of 26 bands, we suggest to use 13 coefficients,
# exactly half of the number of filterbank coefficients.
# if you don't declare ncoeffs, it will be defaulted to half of the number of filterbank coefficients.
# we use a logarithmic rectification.
# the dithered mfcc is suggested, as it's more robust to noise. (taken from AudioFlux Python package)
mfcc = get_mfcc(
    source=melspec,
    ncoeffs = 13,
    rectification = :log, # :log, :cubic_root
    dither = true,
);

# take the fundamental frequency
# we use the NFC method, with a range of 50-400hz.
f0 = get_f0(
    source=stftspec,
    method = :nfc,
	freq_range = (50, 400)
);

# spectral features
# list of spectrals computed:
# spect.centroid
# spect.crest
# spect.entropy
# spect.flatness
# spect.flux
# spect.kurtosis
# spect.rolloff
# spect.skewness
# spect.decrease
# spect.slope
# spect.spread
# also in this case, we use broadband frequency range.
spect = get_spectrals(
    source=stftspec,
    freq_range = (0, round(Int, audio.sr / 2))
);

# collect the features in a single Matrix
# case 1: we want rows as features and columns as frames
audio_features_rows = vcat(
    melspec.spec,
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
    melspec.spec',
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