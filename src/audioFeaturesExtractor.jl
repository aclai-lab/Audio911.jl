"""
Audio911

può essere utilizzata in 2 modi differenti:
#########################################################################################################
1-creare un oggetto audio_obj

audio = audio_features_obj(x, sr)

audio.get_fft()
audio.get_lin_spec()
audio.get_mel_spec()
audio.get_mfcc()
audio.get_spectrals()
audio.get_f0()
audio.get_features(profile, per ora solo :full)

#########################################################################################################
2-utilizzare l'invocazione get_feature per ottenere le features separatamente

costruzione:
get_feature(
    x AbstractFloat file audio mono
    sr Int64 frequenza di campionamento
    feat Symbol audio feature da estrarre
        :fft fast fourier transform (da implementare meglio)
        :lin spettrogramma lineare
        :mel spettrogramma mel
        :mfcc coefficienti mfcc e relative delta e deltadelta
        :spectrals le feature spettrali:
        centroid, crest, decrease, entropy, flatness, flux, kurtosis, rolloff, skewness, slope, spread
        :f0 frequenza fondamentale
    kwargs...

##############################################################################################
paramentri addizionali
sia per l'oggetto audio, che per la chiamata a feature singola

# fft
fft_length::Int64 = 256,
dimensione finestra fft, valori consigliati: 256, 512, 1024

window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
window_length::Int64 = fft_length,
overlap_length::Int64 = round(Int, fft_length * 0.500),
parametri relativi alla finestrazione della fft
di default audio911 usa una finestra di tipo hann, anzichè hamming
con una finestra della stessa dimensione della finestra fft
e un overlap pari alla metà del suo valore

i valori standard sarebbero questi
# window_length::Int64 = Int(round(0.03 * sr)),
# overlap_length::Int64 = Int(round(0.02 * sr)),

window_norm::Bool = false, normalizzazione delle finestre

# spectrum
frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
limiti banda, importantissimi per isolare la porzione di spettro dove si prevede di recuperare l'informazione

spectrum_type::Symbol = :power, # :power, :magnitude
tipo di spettro, di default :power, molto raramente si usa magnitude

# mel
mel_style::Symbol = :htk, # :htk, :slaney
tipo di banco filtro dello spettrogramma mel
la tipologia htk surclassa la tipologia slaney in tutti i nostri esperimenti

mel_bands::Int64 = 26,
numero di bande che compongono lo spettrogramma mel. 26 è il valore di default
ma non c'è un vero e proprio standard di questo valore.
da ricordare che nel caso si voglia utilizzare anche la mfcc,
i coefficienti della mfcc vengono calcolati sulle prime bande dello spettrogramma
quindi una variazione di questo valore comporta anche un diverso funzionamento della mfcc
il cui valore "num_coeff" andrà opportunamente tarato.

filterbank_design_domain::Symbol = :linear,
filterbank_normalization::Symbol = :bandwidth, # :bandwidth, :area, :none
frequency_scale::Symbol = :mel,
l'implementazione corretta di questi paramentri è da completare

# mfcc
num_coeffs::Int64 = 13,
numero delle bande mfcc, vedi sopra

normalization_type::Symbol = :dithered, # :standard, :dithered
paramentro preso da audioflux:
i valori dei coefficienti mfcc, se vanno sotto la soglia di 1e-8
vengono normalizzati a questo valore.
mentre matlab non ha una soglia limite.
dithered > audiofluz, standard > matlab

rectification::Symbol = :log,
log_energy_source::Symbol = :standard, # :standard (after windowing), :mfcc
log_energy_pos::Symbol = :none, #:append, :replace, :none
spesso viene salvato, nella mfcc, il valore del volume in log.
questi parametri definiscono dove viene calcolata e dove salvarla: 
:append viene creato un n-esimo coefficiente, :replace la log energy va a sostituire il primo coefficiente mfcc
se ne sconsiglia comunque l'uso.

delta_window_length::Int64 = 9,
finestra di calcolo della derivata

delta_matrix::Symbol = :transposed, # :standard, :transposed
preso da audioflux che calcola le delta sull'asse delle frequenze anzichè sull'asse temporale
potrebbe sembrare un errore, ma potrebbe anche non esserlo

# spectral
spectral_spectrum::Symbol = :lin, # :lin, :mel
si può scegliere su che spettrogramma calcolare le spectral features: se partendo dal lineare o dal mel

# f0
f0_method::Symbol = :nfc,
f0_range::Tuple{Int64, Int64} = (50, 400),
median_filter_length::Int64 = 1
Questi paramentri sono in fase di studio
"""
################################################################################
#                                audio object                                  #
################################################################################

function audio_features_obj(
        x::AbstractVector{Float64},
        sr::Int64;

        # profile::Symbol = :all,

        # fft
        fft_length::Int64 = 256,
        window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
        window_length::Int64 = fft_length,
        overlap_length::Int64 = round(Int, fft_length * 0.500),
        # window_length::Int64 = Int(round(0.03 * sr)),
        # overlap_length::Int64 = Int(round(0.02 * sr)),
        window_norm::Bool = false,

        # spectrum
        frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
        spectrum_type::Symbol = :power, # :power, :magnitude

        # mel
        mel_style::Symbol = :htk, # :htk, :slaney
        mel_bands::Int64 = 26,
        filterbank_design_domain::Symbol = :linear,
        filterbank_normalization::Symbol = :bandwidth, # :bandwidth, :area, :none
        frequency_scale::Symbol = :mel,

        # mfcc
        num_coeffs::Int64 = 13,
        normalization_type::Symbol = :dithered, # :standard, :dithered
        rectification::Symbol = :log,
        log_energy_source::Symbol = :standard, # :standard (after windowing), :mfcc
        log_energy_pos::Symbol = :none, #:append, :replace, :none
        delta_window_length::Int64 = 9,
        delta_matrix::Symbol = :transposed, # :standard, :transposed

        # spectral
        spectral_spectrum::Symbol = :lin, # :lin, :mel

        # f0
        f0_method::Symbol = :nfc,
        f0_range::Tuple{Int64, Int64} = (50, 400),
        median_filter_length::Int64 = 1
)
    setup = AudioSetup(
        sr = sr,

        # fft
        fft_length = fft_length,
        window_type = window_type,
        window_length = window_length,
        overlap_length = overlap_length,
        window_norm = window_norm,

        # spectrum
        frequency_range = frequency_range,
        spectrum_type = spectrum_type,

        # mel
        mel_style = mel_style,
        mel_bands = mel_bands,
        filterbank_design_domain = filterbank_design_domain,
        filterbank_normalization = filterbank_normalization,
        frequency_scale = frequency_scale,

        # mfcc
        num_coeffs = num_coeffs,
        normalization_type = normalization_type,
        rectification = rectification,
        log_energy_source = log_energy_source,
        log_energy_pos = log_energy_pos,
        delta_window_length = delta_window_length,
        delta_matrix = delta_matrix,

        # spectral
        spectral_spectrum = spectral_spectrum,

        # f0
        f0_method = f0_method,
        f0_range = f0_range,
        median_filter_length = median_filter_length
    )

    # preemphasis
    # zi = 2 * x[1] - x[2]
    # filt!(x, [1.0, -0.97], 1.0, x, [zi])
    # normalize
    # x = x ./ maximum(abs.(x))

    data = AudioData(
        x = x,
    )

    return AudioObj(setup, data)
end

function audio_features_obj(
        x::AbstractVector{T},
        sr::Int64;
        kwargs...
) where {T <: AbstractFloat}
    audio_features_obj(Float64.(x), sr; kwargs...)
end

################################################################################
#                             stand alone functions                            #
################################################################################
function get_fft(setup::AudioSetup, data::AudioData)
    get_fft!(setup, data)

    return data.fft
end

function get_lin_spec(setup::AudioSetup, data::AudioData)
    get_fft!(setup, data)
    lin_spectrogram!(setup, data)

    return data.lin_spectrogram, setup.lin_frequencies
end

function get_mel_spec(setup::AudioSetup, data::AudioData)
    get_fft!(setup, data)
    get_mel_spec!(setup, data)

    return data.mel_spectrogram, setup.mel_frequencies
end

function get_mfcc(setup::AudioSetup, data::AudioData)
    get_fft!(setup, data)
    get_mel_spec!(setup, data)
    get_mfcc!(setup, data)
    get_mfcc_deltas!(setup, data)

    return data.mfcc_coeffs, data.mfcc_delta, data.mfcc_deltadelta
end

function get_spectrals(setup::AudioSetup, data::AudioData)
    get_fft!(setup, data)
    if setup.spectral_spectrum == :lin
        lin_spectrogram!(setup, data)
    elseif setup.spectral_spectrum == :mel
        get_mel_spec!(setup, data)
    else
        error("setup.spectral_spectrum must be :lin or :mel.")
    end
    get_spectrals!(setup, data)

    return [
        data.spectral_centroid,
        data.spectral_crest,
        data.spectral_decrease,
        data.spectral_entropy,
        data.spectral_flatness,
        data.spectral_flux,
        data.spectral_kurtosis,
        data.spectral_rolloff,
        data.spectral_skewness,
        data.spectral_slope,
        data.spectral_spread
    ]
end

function get_f0(setup::AudioSetup, data::AudioData)
    get_f0!(setup, data)

    return data.f0
end

################################################################################
#                        stand alone functions caller                          #
################################################################################
function get_feature(
        x::AbstractVector{Float64},
        sr::Int64,
        feat::Symbol;

        # profile::Symbol = :all,

        # fft
        fft_length::Int64 = 256,
        window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
        window_length::Int64 = fft_length,
        overlap_length::Int64 = round(Int, fft_length * 0.500),
        # window_length::Int64 = Int(round(0.03 * sr)),
        # overlap_length::Int64 = Int(round(0.02 * sr)),
        window_norm::Bool = false,

        # spectrum
        frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
        spectrum_type::Symbol = :power, # :power, :magnitude

        # mel
        mel_style::Symbol = :htk, # :htk, :slaney
        mel_bands::Int64 = 26,
        filterbank_design_domain::Symbol = :linear,
        filterbank_normalization::Symbol = :bandwidth, # :bandwidth, :area, :none
        frequency_scale::Symbol = :mel,

        # mfcc
        num_coeffs::Int64 = 13,
        normalization_type::Symbol = :dithered, # :standard, :dithered
        rectification::Symbol = :log,
        log_energy_source::Symbol = :standard, # :standard (after windowing), :mfcc
        log_energy_pos::Symbol = :none, #:append, :replace, :none
        delta_window_length::Int64 = 9,
        delta_matrix::Symbol = :transposed, # :standard, :transposed

        # spectral
        spectral_spectrum::Symbol = :lin # :lin, :mel
)
    setup = AudioSetup(
        sr = sr,

        # fft
        fft_length = fft_length,
        window_type = window_type,
        window_length = window_length,
        overlap_length = overlap_length,
        window_norm = window_norm,

        # spectrum
        frequency_range = frequency_range,
        spectrum_type = spectrum_type,

        # mel
        mel_style = mel_style,
        mel_bands = mel_bands,
        filterbank_design_domain = filterbank_design_domain,
        filterbank_normalization = filterbank_normalization,
        frequency_scale = frequency_scale,

        # mfcc
        num_coeffs = num_coeffs,
        normalization_type = normalization_type,
        rectification = rectification,
        log_energy_source = log_energy_source,
        log_energy_pos = log_energy_pos,
        delta_window_length = delta_window_length,
        delta_matrix = delta_matrix,

        # spectral
        spectral_spectrum = spectral_spectrum
    )

    data = AudioData(
        x = x,
    )

    calc = Dict([
        :fft => get_fft,
        :lin => get_lin_spec,
        :mel => get_mel_spec,
        :mfcc => get_mfcc,
        :spectrals => get_spectrals,
        :f0 => get_f0])

    return calc[feat](setup, data)
end

function get_feature(
        x::AbstractVector{T},
        sr::Int64,
        feat::Symbol;
        kwargs...
) where {T <: AbstractFloat}
    get_feature(Float64.(x), sr, feat; kwargs...)
end