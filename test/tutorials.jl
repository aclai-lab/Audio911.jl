# ---------------------------------------------------------------------------- #
#                             step-by-step guide                               #
# ---------------------------------------------------------------------------- #

## caricare il pacchetto Audio911
using Audio911

## load an audiofile

# caricare il file audio da cui si vogliono estrarre features.
# il file audio deve essere dei formati .wav, .mp3, .flac o .ogg,
# in caso di file audio stereo, verranno automaticamente convertiti in mono.
# fare sempre molta attenzione alla frequenza di campionamento del file audio:
# se si stanno usando file campionati a 44100, quindi con un range di frequenza 
# da 0 a 22000 Hz, è veramente necessario computare l'estrazione delle features
# per tutta questa banda?
# se poi, successivamente, si imposta un limite di range tipo: `freqrange=[100, 1000]`
# allora è imperativo impostare una corretta frequenza di campionamento che non ecceda
# al range deciso.
# questo porterà ad un notevole aumento delle prestazioni.

# è anche possibile caricare matrici o dataframe di file audio,
# oppure direttamente un file .csv.
# entrambi i casi verranno affrontati nel prossimo tutorial.

# nella cartella `test/test_files` si trovano dei file audio di esempio.
test_file = joinpath(dirname(@__FILE__), "test_files/test.wav")

audio = load(test_file)

# buona pratica sarebbe quella di specificare sempre tutti i parametri:
# nel caso qualcosa andasse storto è più facile risalire al problema.
# altra importantissima buona pratica è quella si specificare sempre il pacchetto da cui 
# proviene la funzione.
# in questo caso è fondamentale, visto che la funzione `load` è definita in tantissimi
# altri pacchetti di julia.
# quindi il metodo consigliato per caricare correttamente un file audio è il seguente:

audio = Audio911.load(test_file; sr=8000, format=Float32, norm=true)

# si richiede che l'audio venga caricato nel modo seguente:
# - se la frequenza di campionamento originale del file è differente, ricampionarlo a 8000 Hz
# - il formato può essere `Float32` o `Float64`.
# - si richiede di normalizzare l'audio (aumentare, se possibile, il volume al massimo)

## from time-domain to frequency domain using fast Fourier transform

# la trasformata di Fourier necessita che l'audio sia dapprima diviso in frames

audioframes = Audio911.Frames(audio; winsize=512, winstep=492, type=hanning)

# quindi è possibile calcolare la stft passandogli i frames

stft_spec = Audio911.Stft(audioframes, nfft=1024, spectrum=power)

stft_spec = Audio911.Stft(audio; nfft=1024, winsize=512, winstep=492, type=hanning, spectrum=power)

# Stft{Float32}
#   Dimensions:
#     Frames:          34
#     Frequency bins:  513
#   Configuration:
#     Sample rate:     8000 Hz
#     FFT size:        1024
#     Window size:     512 samples
#     Overlap:         20 samples
#     Hop size:        492 samples
#     Spectrum type:   power
#   Frequency:
#     Range:           0.0 - 4000.0 Hz

