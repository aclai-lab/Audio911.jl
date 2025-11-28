# ---------------------------------------------------------------------------- #
#                               filterbank setup                               #
# ---------------------------------------------------------------------------- #
struct FBankSetup <: AbstractSetup
    sr        :: Int64
    nbands    :: Int64
    scale     :: Symbol        # :htk, :slaney, :erb, :bark, :semitones
    norm      :: Base.Callable # bandwidth, area, none_norm
    freqrange :: FreqRange
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
        filterbank :: AbstractArray{T},
        filtfreq   :: AbstractVector{T},
        bw         :: AbstractVector{T},
        sr         :: Int64,
        nbands     :: Int64,
        scale      :: Symbol,
        norm       :: Base.Callable,
        freqrange  :: FreqRange,
    ) where {T<:AudioData}
        new{T}(filterbank, filtfreq, bw, FBankSetup(sr, nbands, scale, norm, freqrange))
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    get_data(f::FBank) -> AbstractArray

Get the filter bank, returned as an M-by-N matrix,
where M is the number of bands,
and N is the number of frequency points of a one-sided spectrum.
"""
@inline get_data(f::FBank)          = f.fbank

"""
    get_freq(f::FBank) -> AbstractVector

Get the center frequencies of bandpass filters in Hz.
```
"""
@inline get_freq(f::FBank)          = f.freq

"""
    get_bandwidth(f::FBank) -> AbstractVector

Get the bandwidth of bandpass filters in Hz.
"""
@inline get_bandwidth(f::FBank)     = f.bw

"""
    get_sr(f::FBank) -> Int64

Get the sampling rate used to design the filterbank.
"""
@inline get_sr(f::FBank)            = f.setup.sr

"""
    get_nbands(f::FBank) -> Int64

Get the number of filter bands in the filterbank.
"""
@inline get_nbands(f::FBank)        = f.setup.nbands

"""
    get_scale(f::FBank) -> Symbol

Get the frequency scale used in the filterbank design.
"""
@inline get_scale(f::FBank)         = f.setup.scale

"""
    get_norm(f::FBank) -> Base.Callable

Get the normalization function used in the filterbank.
"""
@inline get_norm(f::FBank)          = f.setup.norm

"""
    get_freqrange(f::FBank) -> FreqRange

Get the frequency range covered by the filterbank.
"""
@inline get_freqrange(f::FBank)     = f.setup.freqrange

function Base.show(io::IO, fb::FBank)
    println(io, "FBank:")
    println(io, "  Sample Rate: $(get_sr(fb)) Hz")
    println(io, "  Number of Bands: $(get_nbands(fb))")
    println(io, "  Scale: $(get_scale(fb))")
    println(io, "  Normalization: $(get_norm(fb))")
    println(io, "  Frequency Range: $(get_low(get_freqrange(fb)))-$(get_hi(get_freqrange(fb))) Hz")
    print(io, "  Filterbank Size: $(size(get_data(fb)))")
end

function Base.show(io::IO, ::MIME"text/plain", fb::FBank)
    show(io, fb)
end

# ---------------------------------------------------------------------------- #
#                         scale convertions functions                          #
# ---------------------------------------------------------------------------- #
const htk(hz::Union{StepRangeLen{T},Vector{T}} where T<:AudioData) = @. 2595 * log10(1 + hz / 700)

const htk(::Type{T}, hz::FreqRange, nbands::Int64) where {T<:AudioData} = begin
    melrange = @. 2595 * log10(1 + T.(hz) / 700)
    melvec = LinRange(get_low(melrange), get_hi(melrange), nbands + 2)  
    return @. 700 * (exp10(melvec / 2595) - 1)
end

const slaney(hz::Union{StepRangeLen{T},Vector{T}} where T<:AudioData) = begin
    lin_step = 200 / 3
    return @. ifelse(hz < 1000, hz / lin_step,
        log(hz * 0.001) / (log(6.4) / 27) + (1000 / lin_step))  
end

const slaney(::Type{T}, hz::FreqRange, nbands::Int64) where {T<:AudioData} = begin
    lin_step = T(200 / 3)
    cp_mel = T(1000 / lin_step)
    hz_T = T.(hz)
    melrange = @. ifelse(hz_T < 1000, hz_T / lin_step,
        log(hz_T * 0.001) / (log(6.4) / 27) + (1000 / lin_step))  
    melvec = LinRange(get_low(melrange), get_hi(melrange), nbands + 2)        
    return @. ifelse(melvec < cp_mel, melvec * lin_step,
        1000 * exp(log(6.4) / 27 * (melvec - cp_mel)))
end

const bark(hz::Union{StepRangeLen{T},Vector{T}} where T<:AudioData) = @. 26.81 * hz / (1960 + hz) - 0.53

const bark(::Type{T}, hz::FreqRange, nbands::Int64) where {T<:AudioData} = begin
    hz_T = T.(hz)
    bark_val = @. 26.81 * hz_T / (1960 + hz_T) - 0.53
    barkrange = map(x -> x < 2 ? T(0.85) * x + T(0.3) : x > T(20.1) ? T(1.22) * x - T(4.422) : x, bark_val)
    barkvec = LinRange(get_low(barkrange), get_hi(barkrange), nbands + 2)
    bark2 = map(x -> x < 2 ? (x - T(0.3)) / T(0.85) : x > T(20.1) ? (x + T(0.22) * T(20.1)) / T(1.22) : x, barkvec)
    return @. 1960 * (bark2 + 0.53) / (26.28 - bark2)
end

# ---------------------------------------------------------------------------- #
#                                 normalization                                #
# ---------------------------------------------------------------------------- #
const area      = (filterbank, bw) -> sum(filterbank, dims=2)
const bandwidth = (filterbank, bw) -> bw / 2
const none_norm = (filterbank, bw) -> 1

function normalize!(
    filterbank :: AbstractArray{T},
    norm_func  :: Function,
    bw         :: AbstractVector{T}
) where {T<:AudioData}
    weight_per_band = norm_func(filterbank, bw)
    @. filterbank  /= (weight_per_band + (weight_per_band == 0))
end

# ---------------------------------------------------------------------------- #
#                           gammatone coefficients                             #
# ---------------------------------------------------------------------------- #
function compute_gammatone_coeffs(sr::Int64, bands::AbstractVector{<:AudioData})
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

    nbands = length(bands)
    b1 = @. -2 * cos_e
    b2 = @. exp(-2 * filt * t)

    cat([[
            t/gain[ind] a11[ind]/gain[ind] 0 1 b1[ind] b2[ind]
            t           a12[ind]           0 1 b1[ind] b2[ind]
            t           a13[ind]           0 1 b1[ind] b2[ind]
            t           a14[ind]           0 1 b1[ind] b2[ind]
        ] for ind in 1:nbands
    ]..., dims=3)
end

# ---------------------------------------------------------------------------- #
#                           design filterbank matrix                           #
# ---------------------------------------------------------------------------- #
function auditory_fbank(
        sr            :: Int64;
        sfreq         :: Union{StepRangeLen{<:AudioData},Nothing}=nothing,
        nfft          :: Int64=512,
        nbands        :: Int64=26,
        scale         :: Function=htk,       # :htk, :slaney, :bark
        norm          :: Function=bandwidth, # area, bandwidth, or none_norm
        domain        :: Symbol=:linear,
        freqrange     :: FreqRange=(0, round(Int, sr / 2))
)
    if isnothing(sfreq)
	    spec_length = get_onesided_stft_range(nfft)[end]
	    sfreq = (0:spec_length - 1) .* (sr / nfft)
    end

    T = eltype(sfreq)
    domain == :warped && (linfq = (0:nfft - 1) .* (sr / nfft))
    band_edges = scale(T, freqrange, nbands)

    filtfreq = @view band_edges[2:(end - 1)]
    nbands = length(filtfreq)

    p = [findfirst(sfreq .> edge) for edge in band_edges]
    isnothing(p[end]) ? p[end] = length(sfreq) : nothing

    # create triangular filters for each band
    filterbank = zeros(T, nbands, length(sfreq))

    # bandwidth
    bw = @views band_edges[3:end] .- band_edges[1:(end - 2)]

    # apply warping transformation if domain is warped
    domain == :warped && (band_edges = scale(band_edges); sfreq = scale(sfreq))

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
    
    FBank(filterbank, filtfreq, bw, sr, nbands, nameof(scale), norm, freqrange)
end

auditory_fbank(s::AbstractSpectrogram; kwargs...) =
    auditory_fbank(get_sr(s); sfreq=get_freq(s), nfft=get_nfft(s), kwargs...)

function gammatone_fbank(
        sr            :: Int64;
        nfft          :: Int64=512,
        nbands        :: Int64=26,
        norm          :: Function=bandwidth, # area, bandwidth, or none_norm
        freqrange     :: FreqRange=(0, round(Int, sr / 2))
)
    erbrange = @. log(10) * 1000 / (24.673 * 4.368) * log10(1 + 0.004368 * freqrange)
    erb      = LinRange(get_low(erbrange), get_hi(erbrange), nbands)
    filtfreq = @. (10 ^ (erb / (log(10) * 1000 / (24.673 * 4.368))) - 1) / 0.004368
    coeffs   = compute_gammatone_coeffs(sr, filtfreq)

    iirfreqz = (b, a, n)-> fft([b; zeros(n - length(b))]) ./ fft([a; zeros(n - length(a))])
    sosfilt  = (c, n)   -> reduce((x, y) -> x .* y, map(row -> iirfreqz(row[1:3], row[4:6], n), eachrow(c)))
    applysos = (i)      -> abs.(sosfilt(coeffs[:, :, i], nfft))

    filterbank = hcat(map(applysos, 1:nbands)...)'
    bw = 1.019 * 24.7 * (0.00437 * filtfreq .+ 1)

    (norm != :none) && normalize!(filterbank, norm, bw)

    @views rem(nfft, 2) == 0 ? filterbank[:, 2:(nfft ÷ 2)] .*= 2 : filterbank[:, 2:(nfft ÷ 2 + 1)] .*= 2
    filterbank = @view filterbank[:, 1:(nfft ÷ 2 + 1)]
    
    FBank(filterbank, filtfreq, bw, sr, nbands, :erb, norm, freqrange)
end

gammatone_fbank(s::AbstractSpectrogram; kwargs...) =
    gammatone_fbank(get_sr(s); nfft=get_nfft(s), kwargs...)
