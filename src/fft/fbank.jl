# ---------------------------------------------------------------------------- #
#                               filterbank setup                               #
# ---------------------------------------------------------------------------- #
struct FBankSetup <: AbstractSetup
    sr            :: Int64
    nbands        :: Int64
    scale         :: Symbol        # :htk, :slaney, :erb, :bark, :semitones
    norm          :: Base.Callable # bandwidth, area, none_norm
    freqrange     :: FreqRange
    semitonerange :: FreqRange
end

# ---------------------------------------------------------------------------- #
#                                filterbank data                               #
# ---------------------------------------------------------------------------- #
struct FBank{T<:AudioData} <: AbstractFBank
    fbank :: AbstractArray{T}
	freq  :: AbstractVector{T}
    bw    :: AbstractVector{T}
    setup :: FBankSetup

    function FBank(
        filterbank    :: AbstractArray{T},
        filt_freq     :: AbstractVector{T},
        bw            :: AbstractVector{T},
        sr            :: Int64,
        nbands        :: Int64,
        scale         :: Symbol,
        norm          :: Base.Callable,
        freqrange     :: FreqRange,
        semitonerange :: FreqRange
    ) where {T<:AudioData}
        new{T}(filterbank, filt_freq, bw, FBankSetup(sr, nbands, scale, norm, freqrange, semitonerange))
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
@inline get_data(f::FBank)          = f.fbank
@inline get_freq(f::FBank)          = f.freq
@inline get_bandwidth(f::FBank)     = f.bw
@inline get_sr(f::FBank)            = f.setup.sr
@inline get_nbands(f::FBank)        = f.setup.nbands
@inline get_scale(f::FBank)         = f.setup.scale
@inline get_norm(f::FBank)          = f.setup.norm
@inline get_freqrange(f::FBank)     = f.setup.freqrange
@inline get_semitonerange(f::FBank) = f.setup.semitonerange

function Base.show(io::IO, fb::FBank)
    println(io, "FBank:")
    println(io, "  Sample Rate: $(get_sr(fb)) Hz")
    println(io, "  Number of Bands: $(get_nbands(fb))")
    println(io, "  Scale: $(get_scale(fb))")
    println(io, "  Normalization: $(get_norm(fb))")
    println(io, "  Frequency Range: $(get_low(get_freqrange(fb)))-$(get_hi(get_freqrange(fb))) Hz")
    # if fb.setup.scale in [:semitones, :tuned_semitones]
    #     println(io, "  Semitone Range: $(get_low(fb.setup.semitonerange))-$(get_hi(fb.setup.semitonerange))")
    # end
    print(io, "  Filterbank Size: $(size(get_data(fb)))")
end

function Base.show(io::IO, ::MIME"text/plain", fb::FBank)
    show(io, fb)
end

# ---------------------------------------------------------------------------- #
#                         scale convertions functions                          #
# ---------------------------------------------------------------------------- #
function hz2mel(
    hz::Union{FreqRange,StepRangeLen{T},Vector{T}},
    style::Symbol
) where {T<:AudioData}
    style == :htk    && return @. 2595 * log10(1 + hz / 700)
    style == :slaney && begin
        lin_step = 200 / 3
        return @. ifelse(hz < 1000, hz / lin_step,
            log(hz * 0.001) / (log(6.4) / 27) + (1000 / lin_step))
    end
    error("Unknown style ($style).")
end

function mel2hz(mel_range::ScaleRange, n_bands::Int64, style::Symbol)
    mel = LinRange(get_low(mel_range), get_hi(mel_range), n_bands + 2)
    style == :htk    && return @. 700 * (exp10(mel / 2595) - 1)
    style == :slaney && begin
        lin_step = 200 / 3
        cp_mel = 1000 / lin_step
        return @. ifelse(
            mel < cp_mel, mel * lin_step, 1000 * exp(log(6.4) / 27 * (mel - cp_mel)))
    end
    error("Unknown style ($style).")
end

function hz2erb(hz::FreqRange)
    @. log(10) * 1000 / (24.673 * 4.368) * log10(1 + 0.004368 * hz)
end

function erb2hz(erb_range::ScaleRange, n_bands::Int64)
    erb = LinRange(get_low(erb_range), get_hi(erb_range), n_bands)
    @. (10 ^ (erb / (log(10) * 1000 / (24.673 * 4.368))) - 1) / 0.004368
end

function hz2bark(hz::FreqRange)
    bark = @. 26.81 * hz / (1960 + hz) - 0.53
    map(x -> x < 2 ? 0.85 * x + 0.3 : x > 20.1 ? 1.22 * x - 4.422 : x, bark)
end

function bark2hz(bark_range::ScaleRange, n_bands::Int64)
    bark = LinRange(get_low(bark_range), get_hi(bark_range), n_bands + 2)
    bark = map(x -> x < 2 ? (x - 0.3) / 0.85 : x > 20.1 ? (x + 0.22 * 20.1) / 1.22 : x, bark)
    @. 1960 * (bark + 0.53) / (26.28 - bark)
end

function hz2semitone(hz::FreqRange)
    get_low(hz) == 0 ? hz = (20, get_hi(hz)) : nothing
    @. 12 * log2(hz)
end

function semitone2hz(st_range::ScaleRange, nbands::Int64)
    st_range_vec = get_low(st_range) .+ 
        collect(0:(nbands+1)) / (nbands+1) * (get_hi(st_range) - get_low(st_range))
    @. 2 ^ (st_range_vec / 12)
end

# ---------------------------------------------------------------------------- #
#                                 normalization                                #
# ---------------------------------------------------------------------------- #
# normalization functions
const area      = (filterbank, bw) -> sum(filterbank, dims=2)
const bandwidth = (filterbank, bw) -> bw / 2
const none_norm = (filterbank, bw) -> 1

function normalize!(filterbank::AbstractArray{Float64}, norm_func::Function, bw::AbstractVector{Float64})
    weight_per_band = norm_func(filterbank, bw)
    @. filterbank /= (weight_per_band + (weight_per_band == 0))
end

# ---------------------------------------------------------------------------- #
#                                   gammatone                                  #
# ---------------------------------------------------------------------------- #
function compute_gammatone_coeffs(sr::Int64, bands::AbstractVector{Float64})
    t = 1 / sr
    erb = @. bands / 9.26449 + 24.7
    filt = 1.019 * 2π * erb

    sqrt_plus = √(3 + 2^(3/2))
    sqrt_minus = √(3 - 2^(3/2))
    img = @. 2im * π * bands * t

    exp_f = @. exp(filt * t)
    exp_im = @. exp(2img)
    exp_it = @. -2 * exp_im * t
    exp_t = @. 2 * exp(-(filt * t) + img) * t

    cos_t = @. cos(2 * bands * π * t)
    sin_t = @. sin(2 * bands * π * t)
    cos_e = @. cos_t / exp_f
    sin_e = @. sin_t / exp_f

    a11 = @. -(t * cos_e + sqrt_plus * t * sin_e)
    a12 = @. -(t * cos_e - sqrt_plus * t * sin_e)
    a13 = @. -(t * cos_e + sqrt_minus * t * sin_e)
    a14 = @. -(t * cos_e - sqrt_minus * t * sin_e)

    gain = @. abs((
        (exp_it + exp_t * (cos_t - sqrt_minus * sin_t)) *
        (exp_it + exp_t * (cos_t + sqrt_minus * sin_t)) *
        (exp_it + exp_t * (cos_t - sqrt_plus * sin_t)) *
        (exp_it + exp_t * (cos_t + sqrt_plus * sin_t))) /
        (-2 / exp(2 * filt * t) - 2 * exp_im + 2 * (1 + exp_im) / exp_f)^4
    )

    n_bands = length(bands)
    b1 = @. -2 * cos_e
    b2 = @. exp(-2 * filt * t)

    cat([
            [
                t/gain[ind] a11[ind]/gain[ind] 0 1 b1[ind] b2[ind]
                t           a12[ind]           0 1 b1[ind] b2[ind]
                t           a13[ind]           0 1 b1[ind] b2[ind]
                t           a14[ind]           0 1 b1[ind] b2[ind]
            ] for ind in 1:n_bands
        ]..., dims=3)
end

# ---------------------------------------------------------------------------- #
#                           design filterbank matrix                           #
# ---------------------------------------------------------------------------- #
function FBank(
        sr            :: Int64;
        sfreq         :: Union{StepRangeLen{<:Float64},Nothing}=nothing,
        nfft          :: Int64=512,
        nbands        :: Int64=26,
        scale         :: Symbol=:htk, # :htk, :slaney, :erb, :bark, :semitones
        norm          :: Function=bandwidth, # area, bandwidth, or none_norm
        domain        :: Symbol=:linear,
        freqrange     :: Tuple{Int64, Int64}=(0, round(Int, sr / 2)),
        semitonerange :: Tuple{Int64, Int64}=(200, 700),
)
    if isnothing(sfreq)
	    spec_length = get_onesided_stft_range(nfft)[end]
	    sfreq = (0:spec_length - 1) .* (sr / nfft)
    end

    domain == :warped && (linfq = (0:nfft - 1) .* (sr / nfft))

    if scale == :erb
        erb_range = hz2erb(freqrange)
        filt_freq = erb2hz(erb_range, nbands)
        coeffs    = compute_gammatone_coeffs(sr, filt_freq)

        iirfreqz  = (b, a, n)-> fft([b; zeros(n - length(b))]) ./ fft([a; zeros(n - length(a))])
        sosfilt   = (c, n)   -> reduce((x, y) -> x .* y, map(row -> iirfreqz(row[1:3], row[4:6], n), eachrow(c)))
        apply_sos = (i)      -> abs.(sosfilt(coeffs[:, :, i], nfft))

        filterbank = hcat(map(apply_sos, 1:nbands)...)'
        bw = 1.019 * 24.7 * (0.00437 * filt_freq .+ 1)

        (norm != :none) && normalize!(filterbank, norm, bw)

        @views rem(nfft, 2) == 0 ? filterbank[:, 2:(nfft ÷ 2)] .*= 2 : filterbank[:, 2:(nfft ÷ 2 + 1)] .*= 2
        filterbank = @view filterbank[:, 1:(nfft ÷ 2 + 1)]

    else
        if (scale == :htk || scale == :slaney)
            mel_range  = hz2mel(freqrange, scale)
            band_edges = mel2hz(mel_range, nbands, scale)
        elseif  scale == :bark
            bark_range = hz2bark(freqrange)
            band_edges = bark2hz(bark_range, nbands)
        elseif scale == :semitones
            st_range   = hz2semitone(freqrange)
            band_edges = semitone2hz(st_range, nbands)
        else
            error("Unknown filterbank frequency scale '$scale', available scales are: :htk, :slaney, :erb, :bark, :semitones, , :semitones_tuned.")
        end

        filt_freq = @view band_edges[2:(end - 1)]
        nbands = length(filt_freq)

        p = [findfirst(sfreq .> edge) for edge in band_edges]
        isnothing(p[end]) ? p[end] = length(sfreq) : nothing

        # create triangular filters for each band
        filterbank = zeros(nbands, length(sfreq))

        # bandwidth
        bw = @views band_edges[3:end] .- band_edges[1:(end - 2)]

        # Apply warping transformation if domain is warped
        if domain == :warped
            band_edges, sfreq = if scale == :htk || scale == :slaney
                hz2mel(band_edges, scale), hz2mel(sfreq, scale)
            elseif scale == :bark
                hz2bark(band_edges), hz2bark(sfreq)
            else
                band_edges, sfreq
            end
        end

        for k in 1:nbands
            # rising side of triangle
            rise_range = p[k]:(p[k + 1] - 1)
            @views @. filterbank[k, rise_range] =
                (sfreq[rise_range] - band_edges[k]) / (band_edges[k + 1] - band_edges[k])
            # falling side of triangle
            fall_range = p[k + 1]:(p[k + 2] - 1)
            @views @. filterbank[k, fall_range] =
                (band_edges[k + 2] - sfreq[fall_range]) / (band_edges[k + 2] - band_edges[k + 1])
        end

        # normalization
        (norm != :none) && normalize!(filterbank, norm, bw)
    end
    
    FBank(filterbank, filt_freq, bw, sr, nbands, scale, norm, freqrange, semitonerange)
end

FBank(s::AbstractSpectrogram; kwargs...) =
    FBank(get_sr(s); sfreq=get_freq(s), nfft=get_nfft(s), kwargs...)
