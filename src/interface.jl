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

# ---------------------------------------------------------------------------- #
#                                  Audio911                                    #
# ---------------------------------------------------------------------------- #
function afe(x::Union{String, AbstractVector{Float64}}; kwargs...)
    aparams = (AUDIO_PARAMS[k] => v for (k, v) in kwargs if k in keys(AUDIO_PARAMS))
    if isa(x, String)
        audio = load_audio(; source=x, aparams...)
    else
        audio =Audio(x, kwargs[:sr])
    end

    featset = get(kwargs, :featset, ()) |> x -> isa(x, Symbol) ? (x,) : x

    stftspec = get_stft(audio; (STFT_PARAMS[k] => v for (k, v) in kwargs if k in keys(STFT_PARAMS))...)

    :lin in featset && begin
        linspec = get_linspec(; source=stftspec, (LIN_PARAMS[k] => v for (k, v) in kwargs if k in keys(LIN_PARAMS))...)
    end

    melfb = get_melfb(; source=stftspec, (MEL_PARAMS[k] => v for (k, v) in kwargs if k in keys(MEL_PARAMS))...)

    :get_only_freqs in featset && return melfb.data.freq

    melspec = get_melspec(; source=stftspec, fbank=melfb, (MELSPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(MELSPEC_PARAMS))...)

    :mfcc in featset && begin
        mfcc = get_mfcc(; source=melspec, (MFCC_PARAMS[k] => v for (k, v) in kwargs if k in keys(MFCC_PARAMS))...)
        :deltas in featset && begin
            deltas = get_deltas(; source=mfcc, (DELTA_PARAMS[k] => v for (k, v) in kwargs if k in keys(DELTA_PARAMS))...)
        end
    end

    :f0 in featset && begin
        f0 = get_f0(; source=stftspec, (F0_PARAMS[k] => v for (k, v) in kwargs if k in keys(F0_PARAMS))...)
    end

    spect = get_spectrals(; source=stftspec, (SPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(SPEC_PARAMS))...)

    :cwt in featset && begin
        cwtfb = get_cwtfb(; source=audio, (CWT_PARAMS[k] => v for (k, v) in kwargs if k in keys(CWT_PARAMS))...)
        cwtspec = get_cwt(; source=audio, fbank=cwtfb, (CWTSPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(CWTSPEC_PARAMS))...)
    end

    hcat(
        filter(!isnothing, [
            :lin in featset ? linspec.data.spec' : nothing,
            melspec.data.spec',
            :mfcc in featset ? mfcc.data.spec' : nothing,
            :deltas in featset ? hcat(deltas.data.dspec', deltas.data.ddspec') : nothing,
            :f0 in featset ? f0.data.f0 : nothing,
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
        ])...
    )  
end

function mfcc(audio::Audio; kwargs...)
    stftspec = get_stft(; source=audio, (STFT_PARAMS[k] => v for (k, v) in kwargs if k in keys(STFT_PARAMS))...)
    melfb = get_melfb(; source=stftspec, (MEL_PARAMS[k] => v for (k, v) in kwargs if k in keys(MEL_PARAMS))...)
    melspec = get_melspec(; source=stftspec, fbank=melfb, (MELSPEC_PARAMS[k] => v for (k, v) in kwargs if k in keys(MELSPEC_PARAMS))...)
    mfcc = get_mfcc(; source=melspec, (MFCC_PARAMS[k] => v for (k, v) in kwargs if k in keys(MFCC_PARAMS))...)

    mfcc.data.spec'
end