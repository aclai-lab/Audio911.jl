# ---------------------------------------------------------------------------- #
#                                backup functs                                 #
# ---------------------------------------------------------------------------- #
# function cep2spec(cep::Matrix{Float64}, nfreq::Int=21, type::Int=2)
#     ncep, ncol = size(cep)

#     # Make the DCT matrix
#     dctm = zeros(ncep, nfreq)
#     idctm = zeros(nfreq, ncep)

#     if type == 2 || type == 3
#         # This is the orthogonal one, so inv matrix is same as fwd matrix
#         for i in 1:ncep
#             dctm[i,:] = cos.((i-1) * (1:2:(2*nfreq-1)) / (2*nfreq) * π) * sqrt(2/nfreq)
#         end
        
#         if type == 2
#             # Make it unitary! (but not for HTK type 3)
#             dctm[1,:] ./= sqrt(2)
#         else
#             dctm[1,:] ./= 2
#         end
#         idctm = dctm'
#     elseif type == 4  # Type 1 with implicit repetition of first, last bins
#         # Reconstruct the middle nfreq rows of an nfreq+2 row idctm
#         for i in 1:ncep
#             # 2x to compensate for fact that only getting +ve freq half
#             idctm[:,i] = 2 * cos.((i-1) * (1:nfreq) / (nfreq+1) * π)
#         end
#         # Fixup 'non-repeated' basis fns
#         idctm[:, [1, ncep]] ./= 2
#     else  # dpwe type 1 - idft of cosine terms
#         for i in 1:ncep
#             # 2x to compensate for fact that only getting +ve freq half
#             idctm[:,i] = 2 * cos.((i-1) * (0:(nfreq-1)) / (nfreq-1) * π)
#         end
#         # Fixup 'non-repeated' basis fns
#         idctm[:, [1, ncep]] .*= 0.5
#     end

#     spec = exp.(idctm * cep)
    
#     return spec, idctm
# end

function hz2mel(hz::Tuple{Int64, Int64}, style::Symbol)
    style == :mel_htk && return @. 2595 * log10(1 + hz / 700)
    style == :mel_slaney && begin
        lin_step = 200 / 3
        return @. ifelse(hz < 1000, hz / lin_step,
            log(hz * 0.001) / (log(6.4) / 27) + (1000 / lin_step))
    end
    error("Unknown style ($style).")
end

function mel2hz(mel_range::Tuple{Float64, Float64}, n_bands::Int64, style::Symbol)
    mel = LinRange(mel_range[1], mel_range[2], n_bands + 2)
    style == :mel_htk && return @. 700 * (exp10(mel / 2595) - 1)
    style == :mel_slaney && begin
        lin_step = 200 / 3
        cp_mel = 1000 / lin_step
        return @. ifelse(
            mel < cp_mel, mel * lin_step, 1000 * exp(log(6.4) / 27 * (mel - cp_mel)))
    end
    error("Unknown style ($style).")
end

function mfcc2mel(mfcc::AbstractArray{<:AbstractFloat}, nbands::Int64, ncoeffs::Int64)
    # Make the DCT matrix
    dctm = zeros(ncoeffs, nbands)
    idctm = zeros(nbands, ncoeffs)

    for i in 1:ncoeffs
        dctm[i,:] = cos.((i-1) * (1:2:(2*nbands-1)) / (2*nbands) * π) * sqrt(2/nbands)
    end
    
    dctm[1,:] ./= sqrt(2)
    idctm = dctm'

    spec = exp.(idctm * mfcc)
    
    return spec, idctm
end

function htk2stft(sr::Int64, nbands::Int64, nfft::Int, freq_range::Tuple{Int64, Int64}, scale::Symbol)
    wts = zeros(nbands, nfft)
    # Center freqs of each FFT bin
    fftfrqs = (0:(nfft÷2)) / nfft * sr

    mel_range = hz2mel(freq_range, scale)
    band_edges = mel2hz(mel_range, nbands, scale)

    binbin = round.(Int, band_edges / sr * (nfft-1))

    for i in 1:nbands
        fs = band_edges[i .+ (0:2)]
        fs = fs[2] * (fs .- fs[2])
        # lower and upper slopes for all bins
        loslope = (fftfrqs .- fs[1]) / (fs[2] - fs[1])
        hislope = (fs[3] .- fftfrqs) / (fs[3] - fs[2])
        # .. then intersect them with each other and zero
        wts[i, 1 .+ (0:(nfft÷2))] = max.(0, min.(loslope, hislope))
    end

    # if !constamp
    # Slaney-style mel is scaled to be approx constant E per channel
    wts = Diagonal(2 ./ (band_edges[2 .+ (1:nbands)] .- band_edges[1:nbands])) * wts
    # end

    # Make sure 2nd half of FFT is zero
    wts[:, (nfft÷2+2):nfft] .= 0

    return wts, band_edges
end

function mel2stft(spec::AbstractArray{<:AbstractFloat}, sr::Int64, nbands::Int64, nfft::Int64, scale::Symbol, freq_range::Tuple{Int64, Int64}, melfb_norm::Symbol)
    if scale == :mel_htk
        wts, _ = htk2stft(sr, nbands, nfft, freq_range, scale)
    # elseif fbtype == :mel_slaney
    #     wts = slan2stft(sr, nbands, nfft, freq_range)
    # elseif fbtype == :bark
    #     wts = bark2stft(sr, nbands, nfft, freq_range, true, true)
    # elseif fbtype == "fcmel"
    #     wts = fft2melmx(nfft, sr, nbands, bwidth, minfreq, maxfreq, true)
    else
        error("scale $scale not recognized")
    end

    # Cut off 2nd half
    wts = wts[:, 1:((nfft÷2)+1)]
    # Cut off 2nd half
    wts = wts[:, 1:((nfft÷2)+1)]
    # Just transpose, fix up
    ww = wts' * wts
    iwts = wts' ./ repeat(max.(StatsBase.mean(diag(ww))/100, sum(ww, dims=2)), 1, nbands)

    # Apply weights
    if melfb_norm == :power
        spec = iwts * spec
    else
        spec = (iwts * sqrt.(spec)).^2
    end

    return spec, wts, iwts
end

function invmfcc(
    mfcc::AbstractArray{<:AbstractFloat};
    sr::Int64,
    nfft::Int64,
    # # win_type = (:hann, :periodic)
    win_length::Int64,
    overlap_length::Int64,
    # stft_norm::Symbol,           # :power, :magnitude, :pow2mag
    nbands::Int64,
    scale::Symbol,                 # :mel_htk, :mel_slaney, :erb, :bark
    melfb_norm::Symbol,                  # :bandwidth, :area, :none
    freq_range::Tuple{Int64, Int64},
    ncoeffs::Int64,
)
    # undo the dct
    inv_mel, _ = mfcc2mel(mfcc, nbands, ncoeffs)
    # undo the auditory spectrum
    inv_stft, _, _ = mel2stft(inv_mel, sr, nbands, nfft, scale, freq_range, melfb_norm)
    # back to waveform
    x = invpowspec(inv_stft, sr, win_length, overlap_length)
end

# ---------------------------------------------------------------------------- #
#                                    debug                                    #
# ---------------------------------------------------------------------------- #
# debug
using Pkg
Pkg.activate("/home/paso/Documents/Aclai/audio-rules2024")
using WAV, Audio911
using DSP, FFTW
using LinearAlgebra, StatsBase

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

sr = 8000
audio = load_audio(file=wavfile, sr=sr, norm=true);

nfft = 512
win_type = (:hann, :periodic)
win_length = 512
overlap_length = 256
stft_norm = :power              # :power, :magnitude, :pow2mag
nbands = 40
scale = :mel_htk                 # :mel_htk, :mel_slaney, :erb, :bark
melfb_norm = :bandwidth                # :bandwidth, :area, :none
freq_range = (0, round(Int, audio.sr / 2))
db_scale = false
ncoeffs = 13
rectification = :log
dither = false

function afe_mel(audio::Audio;
    # -------------------------------- parameters -------------------------------- #
    # audio module
    sr = audio.sr,
    norm = true,
    speech_detection = false,
    # stft module
    nfft = nfft,
    win_type = win_type,
    win_length = win_length,
    overlap_length = overlap_length,
    stft_norm = stft_norm,                      # :power, :magnitude, :pow2mag
    # mel filterbank module
    nbands = nbands,
    scale = scale,                      # :mel_htk, :mel_slaney, :erb, :bark
    melfb_norm = melfb_norm,                 # :bandwidth, :area, :none
    freq_range = freq_range,
    # mel spectrogram module
    db_scale = db_scale,
)
    # --------------------------------- functions -------------------------------- #
    stftspec = get_stft(
        audio=audio,
        nfft=nfft,
        win_type=win_type,
        win_length=win_length,
        overlap_length=overlap_length,
        norm=stft_norm
    );

    # mel filterbank module
    melfb = get_melfb(
        stft=stftspec,
        nbands=nbands,
        scale=scale,
        norm=melfb_norm,
        freq_range=freq_range
    );

    # mel spectrogram module
    get_melspec(
        stft=stftspec,
        fbank=melfb,
        db_scale=db_scale
    )
end

function afe_mfcc(melspec::MelSpec;
    # -------------------------------- parameters -------------------------------- #
    # mfcc module
    ncoeffs = ncoeffs,
    rectification = rectification,                    # :log, :cubic_root
    dither = dither,
)

    # --------------------------------- functions -------------------------------- #
    # mfcc module
    get_mfcc(
        source=melspec,
        ncoeffs=ncoeffs,
        rectification=rectification,
        dither=dither,
    )
end

n_mel_bands = 26
melspec = afe_mel(audio; nbands=n_mel_bands, db_scale=false)
mfccspec = afe_mfcc(melspec)

mfcc = mfccspec.mfcc # debug


########################à
using DSP, FFTW, Random

function invpowspec(inv_stft::AbstractArray{<:AbstractFloat}, sr::Int64, nfft::Int64, win_length::Int64, overlap_length::Int64)
    
    nrow, ncol = size(inv_stft)

    window = DSP.hamming(win_length)
    x = ifft(inv_stft, win_length, overlap_length; nfft=nfft, fs=sr, window=window)

    return real(x)
end

function istft(X::Matrix{Complex{Float64}}, winlen::Int, noverlap::Int; 
    nfft::Int=size(X,1)*2-2, fs::Int=1, window::Vector{Float64}=ones(winlen))

hop = winlen - noverlap
time_len = size(X,2) * hop + winlen
x = zeros(time_len)

for i in 1:size(X,2)
start = (i-1)*hop + 1
ending = start + winlen - 1
x[start:ending] .+= real.(ifft(vcat(X[:,i], conj.(X[end-1:-1:2,i]))))
end

return x
end

