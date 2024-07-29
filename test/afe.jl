using Audio911

# ---------------------------------------------------------------------------- #
#                                    utils                                     #
# ---------------------------------------------------------------------------- #
nan_replacer!(x::AbstractArray{Float64}) = replace!(x, NaN => 0.0)

# ---------------------------------------------------------------------------- #
#                       audio911 audio features extractor                      #
# ---------------------------------------------------------------------------- #
function audio911_extractor(
    # audio module
    wavfile::Union{String, AbstractVector{Float64}};
    sr::Int64=8000,
    norm::Bool=true,
    speech_detection::Bool=false,
    # stft module
    stft_length::Union{Int64, Nothing}=nothing,
    win_type::Tuple{Symbol, Symbol}=(:hann, :periodic),
    win_length::Union{Int64, Nothing}=nothing,
    overlap_length::Union{Int64, Nothing}=nothing,
    stft_norm::Symbol=:power,               # :power, :magnitude, :pow2mag
    # mel filterbank module
    nbands::Int64=26,
    scale::Symbol=:bark,                 # :mel_htk, :mel_slaney, :erb, :bark
    melfb_norm::Symbol=:bandwidth,          # :bandwidth, :area, :none
    freq_range::Union{Tuple{Int64, Int64}, Nothing}=nothing,
    # mel spectrogram module
    db_scale::Bool=false,
    # mfcc module
    ncoeffs::Int64=13,
    rectification::Symbol=:log,             # :log, :cubic_root
    dither::Bool=true,
    # # deltas module
    # d_length = 9,
    # d_matrix = :transposed,                 # :standard, :transposed
    # f0 module
    method::Symbol=:nfc,
    f0_range::Tuple{Int64, Int64}=(50, 400),
    # spectral features module
    spect_range::Union{Tuple{Int64, Int64}, Nothing}=nothing,
)
    # audio module
    audio = load_audio(
        file=wavfile, 
        sr=sr, 
        norm=norm,
    );
    if speech_detection
        audio = speech_detector(audio=audio);
    end

    # stft module
    if isnothing(stft_length)
        stft_length	= audio.sr <= 8000 ? 256 : 512
    end
    if isnothing(win_length)
        win_length = stft_length
    end
    if isnothing(overlap_length)
        overlap_length = round(Int, stft_length / 2)
    end
    stftspec = get_stft(
        audio=audio, 
        stft_length=stft_length,
        win_type=win_type,
        win_length=win_length,
        overlap_length=overlap_length,
        norm=stft_norm
    );

    # mel filterbank module
    if isnothing(freq_range)
        freq_range = (0, round(Int, audio.sr / 2))
    end
    melfb = get_melfb(
        stft=stftspec,
        nbands=nbands,
        scale=scale,
        norm=melfb_norm,
        freq_range=freq_range
    );

    # mel spectrogram module
    melspec =  get_melspec(
        stft=stftspec, 
        fbank=melfb,
        db_scale=db_scale
    );

    # mfcc module
    mfcc = get_mfcc(
        source=melspec,
        ncoeffs=ncoeffs,
        rectification=rectification,
        dither=dither,
    );

    # # deltas module
    # deltas = get_deltas(
    #     source=mfcc,
    #     d_length=d_length,
    #     d_matrix=d_matrix
    # );

    # f0 module
    f0 = get_f0(
        source=stftspec,
        method=method,
        freq_range=f0_range
    );

    # spectral features module
    if isnothing(spect_range)
        spect_range = freq_range
    end
    spect = get_spectrals(
        source=stftspec,
        freq_range=spect_range
    );

    return hcat(
        melspec.spec',
        mfcc.mfcc',
        # deltas.delta',
        # deltas.ddelta',
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
end
