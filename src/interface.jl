# ---------------------------------------------------------------------------- #
#                                 parameters                                   #
# ---------------------------------------------------------------------------- #
const AUDIO_PARAMS = Dict(:sr => :sr, :audio_norm => :norm)
const STFT_PARAMS = Dict(:nfft => :nfft, :fft_wintype => :wintype, :fft_nwin => :nwin, :fft_noverlap => :noverlap, :fft_norm => :norm)
const LIN_PARAMS = Dict(:lin_freqrange => :freqrange, :lin_dbscale => :dbscale)
const MEL_PARAMS = Dict(:mel_nbands => :nbands, :mel_scale => :scale, :mel_norm => :norm, :mel_freqrange => :freqrange, :mel_semitonerange => :semitonerange)
const MELSPEC_PARAMS = Dict(:mel_dbscale => :dbscale)
const MFCC_PARAMS = Dict(:mfcc_ncoeffs => :ncoeffs, :mfcc_rect => :rect, :mfcc_dither => :dither)
const DELTA_PARAMS = Dict(:delta_length => :dlength, :delta_transpose => :transpose)
const F0_PARAMS = Dict(:f0_method => :method, :f0_freqrange => :freqrange, :f0_mflength => :mflength)
const SPEC_PARAMS = Dict(:spec_freqrange => :freqrange)
const CWT_PARAMS = Dict(:cwt_wavelet => :wavelet, :cwt_morseparams => :morseparams, :cwt_vpo => :vpo, :cwt_freqrange => :freqrange)
const CWTSPEC_PARAMS = Dict(:cwt_norm => :norm, :cwt_dbscale => :dbscale)

params(x::Dict, kwargs) = (x[k] => v for (k, v) in kwargs if k in keys(x))

# ---------------------------------------------------------------------------- #
#                                  Audio911                                    #
# ---------------------------------------------------------------------------- #
function audio_features(audio::Audio; kwargs...)
    featset = get(kwargs, :featset, ()) |> x -> isa(x, Symbol) ? (x,) : x

    # stftspec = get_stft(audio; (STFT_PARAMS[k] => v for (k, v) in kwargs if k in keys(STFT_PARAMS))...)
    stftspec = get_stft(audio; params(STFT_PARAMS, kwargs)...)

    :lin in featset && begin
        linspec = get_linspec(stftspec; (LIN_PARAMS[k] => v for (k, v) in kwargs if k in keys(LIN_PARAMS))...)
    end

    (:mel in featset || :mfcc in featset || :get_only_freqs in featset) && begin
        melfb = get_melfb(stftspec; (MEL_PARAMS[k] => v for (k, v) in kwargs if k in keys(MEL_PARAMS))...)

        :get_only_freqs in featset && return melfb.data.freq

        melspec = get_melspec(stftspec; fbank=melfb, (MELSPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(MELSPEC_PARAMS))...)
    end

    :mfcc in featset && begin
        mfcc = get_mfcc(melspec; (MFCC_PARAMS[k] => v for (k, v) in kwargs if k in keys(MFCC_PARAMS))...)
        :deltas in featset && begin
            deltas = get_deltas(mfcc; (DELTA_PARAMS[k] => v for (k, v) in kwargs if k in keys(DELTA_PARAMS))...)
        end
    end

    :f0 in featset && begin
        f0 = get_f0(stftspec; (F0_PARAMS[k] => v for (k, v) in kwargs if k in keys(F0_PARAMS))...)
    end

    :spectrals in featset && begin
        spect = get_spectrals(stftspec; (SPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(SPEC_PARAMS))...)
    end

    :cwt in featset && begin
        cwtfb = get_cwtfb(; source=audio, (CWT_PARAMS[k] => v for (k, v) in kwargs if k in keys(CWT_PARAMS))...)
        cwtspec = get_cwt(; source=audio, fbank=cwtfb, (CWTSPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(CWTSPEC_PARAMS))...)
    end 

    hcat(
    filter(!isnothing, [
        (:lin in featset ? linspec.data.spec' : nothing),
        (:mel in featset ? melspec.data.spec' : nothing),
        (:mfcc in featset ? mfcc.data.spec' : nothing),
        (:deltas in featset ? hcat(deltas.data.dspec', deltas.data.ddspec') : nothing),
        (:f0 in featset ? f0.data.f0 : nothing),
        (:spectrals in featset ? hcat(
            spect.data.centroid,
            spect.data.crest,
            spect.data.entropy,
            spect.data.flatness,
            spect.data.flux,
            spect.data.kurtosis,
            spect.data.rolloff,
            spect.data.skewness,
            spect.data.decrease,
            spect.data.slope,
            spect.data.spread
        ) : nothing)
    ])...
)

end

audio_features(x::AbstractVector{<:AbstractFloat}, sr::Int; kwargs...) = audio_features(Audio(x, sr); kwargs...)
audio_features(x::String, sr::Int; norm::Bool=false, kwargs...) = audio_features(load_audio(x, sr; norm=norm); kwargs...)
