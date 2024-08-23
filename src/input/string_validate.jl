debug = "stft(fft_length=1024), melfb(type=mel), cwt(source=mel_fb), mfcc(source=cwt, ncoeffs=16)"

# --- result ----------------------------------------------------------------- #
# result = audio911("stft(fft_length=1024), melfb(type=erb), cwt(source=mel_fb), mfcc(source=cwt, ncoeffs=16)")

# string >>>
# feats2extract[
#     get_stft(fft_length=1024),
    # melfb pipeline: no extra parameters asked (eg:fft_length), so take the previous fft_length
    get_stft(fft_length=1024),
    get_melfb(type=erb),
    # cwt

"""
i valori vanno specificati nella propria funzione della pipeline.
Se vengono specificati in una funzione successiva (ad esempio: mfcc(fft_length=1024, melbands=40))
allora non devono essere specificati, se incontrate, nelle funzioni precedenti. altrimenti c'è un conflitto.
ci deve essere tassativa coerenza tra i dati.
se uno chiede stft e cwt insieme, e la stft non è della lunghezza precisa della cwt, 
"""