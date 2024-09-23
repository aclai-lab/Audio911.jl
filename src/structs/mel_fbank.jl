# ---------------------------------------------------------------------------- #
#                                mel filterbank                                #
# ---------------------------------------------------------------------------- #
struct MelFbSetup
    nbands::Int
    scale::Symbol # :mel_htk, :mel_slaney, :erb, :bark, :semitones, :tuned_semitones
    norm::Symbol # :bandwidth, :area, :none
    freqrange::Tuple{Int, Int}
    semitonerange::Tuple{Int, Int}
end

struct MelFbData
    fbank::AbstractArray{<:AbstractFloat}
	freq::AbstractVector{<:AbstractFloat}
end

struct MelFb
    sr::Int
    setup::MelFbSetup
    data::MelFbData
end

# ---------------------------------------------------------------------------- #
#                         scale convertions functions                          #
# ---------------------------------------------------------------------------- #
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

function hz2erb(hz::Tuple{Int64, Int64})
    @. log(10) * 1000 / (24.673 * 4.368) * log10(1 + 0.004368 * hz)
end

function erb2hz(erb_range::Tuple{Float64, Float64}, n_bands::Int64)
    erb = LinRange(erb_range[1], erb_range[2], n_bands)
    @. (10 ^ (erb / (log(10) * 1000 / (24.673 * 4.368))) - 1) / 0.004368
end

function hz2bark(hz::Tuple{Int64, Int64})
    bark = @. 26.81 * hz / (1960 + hz) - 0.53
    map(x -> x < 2 ? 0.85 * x + 0.3 : x > 20.1 ? 1.22 * x - 4.422 : x, bark)
end

function bark2hz(bark_range::Tuple{Float64, Float64}, n_bands::Int64)
    bark = LinRange(bark_range[1], bark_range[2], n_bands + 2)
    bark = map(x -> x < 2 ? (x - 0.3) / 0.85 : x > 20.1 ? (x + 0.22 * 20.1) / 1.22 : x, bark)
    @. 1960 * (bark + 0.53) / (26.28 - bark)
end

function hz2semitone(hz::Tuple{Int64, Int64})
    hz[1] == 0 ? hz = (20, hz[2]) : nothing
    @. 12 * log2(hz)
end

function semitone2hz(st_range::Tuple{Float64, Float64}, nbands::Int64)
    st_range_vec = st_range[1] .+ collect(0:(nbands+1)) / (nbands+1) * (st_range[2] - st_range[1])
    @. 2 ^ (st_range_vec / 12)
end

# ---------------------------------------------------------------------------- #
#                                 normalization                                #
# ---------------------------------------------------------------------------- #
function normalize!(filterbank::AbstractArray{Float64}, norm::Symbol, bw::AbstractVector{Float64})
    norm_funcs = Dict(
        :area => () -> sum(filterbank, dims = 2),
        :bandwidth => () -> bw / 2,
        :none => () -> 1
    )

    weight_per_band = get(norm_funcs, norm, norm_funcs[:none])()
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

# -------------------------------------------------------------------------- #
#                                  f0 log mel                                #
# -------------------------------------------------------------------------- #
function calc_f0(source::AbstractVector{<:AbstractFloat}, sfreq::StepRangeLen{<:AbstractFloat}, f0_semitone_range::Tuple{Int, Int},)
	x1 = findfirst((x)-> x >= f0_semitone_range[1], sfreq)
	x2 = findfirst((x)-> x >= f0_semitone_range[2], sfreq)

	peak_pos = argmax(source[x1:x2-1])
	sfreq[x1+peak_pos-1]
end

# ---------------------------------------------------------------------------- #
#                           design filterbank matrix                           #
# ---------------------------------------------------------------------------- #
function _get_melfb(;
        sr::Int,
        source::AbstractVector{<:AbstractFloat},
        sfreq::StepRangeLen{<:AbstractFloat},
        nfft::Int,
        nbands::Int = 26,
        scale::Symbol = :mel_htk, # :mel_htk, :mel_slaney, :erb, :bark, :semitones, :tuned_semitones
        norm::Symbol = :bandwidth,  # :bandwidth, :area, :none
        freqrange::Tuple{Int, Int} = (0, round(Int, sr / 2)),
        semitonerange::Tuple{Int, Int} = (200, 700),
)
    if scale == :erb
        erb_range = hz2erb(freqrange)
        filter_freq = erb2hz(erb_range, nbands)
        coeffs = compute_gammatone_coeffs(sr, filter_freq)

        iirfreqz = (b, a, n) -> fft([b; zeros(n - length(b))]) ./ fft([a; zeros(n - length(a))])
        sosfilt = (c, n) -> reduce((x, y) -> x .* y, map(row -> iirfreqz(row[1:3], row[4:6], n), eachrow(c)))
        apply_sosfilt = (i) -> abs.(sosfilt(coeffs[:, :, i], nfft))

        filterbank = hcat(map(apply_sosfilt, 1:nbands)...)'
        bw = 1.019 * 24.7 * (0.00437 * filter_freq .+ 1)

        (norm != :none) && normalize!(filterbank, norm, bw)

        rem(nfft, 2) == 0 ? filterbank[:, 2:(nfft ÷ 2)] .*= 2 : filterbank[:, 2:(nfft ÷ 2 + 1)] .*= 2
        filterbank = filterbank[:, 1:(nfft ÷ 2 + 1)]

    else
        if (scale == :mel_htk || scale == :mel_slaney)
            mel_range = hz2mel(freqrange, scale)
            band_edges = mel2hz(mel_range, nbands, scale)
        elseif  scale == :bark
            bark_range = hz2bark(freqrange)
            band_edges = bark2hz(bark_range, nbands)
        elseif scale == :semitones
            st_range = hz2semitone(freqrange)
            band_edges = semitone2hz(st_range, nbands)
        elseif scale == :tuned_semitones
            freqrange = (round(Int, calc_f0(source, sfreq, semitonerange)), freqrange[2])
            st_range = hz2semitone(freqrange)
            band_edges = semitone2hz(st_range, nbands)
        else
            error("Unknown filterbank frequency scale '($scale)', available scales are: :mel_htk, :mel_slaney, :erb, :bark, :semitones, , :semitones_tuned.")
        end

        filter_freq = band_edges[2:(end - 1)]
        nbands = length(filter_freq)

        p = [findfirst(sfreq .> edge) for edge in band_edges]
        isnothing(p[end]) ? p[end] = length(sfreq) : nothing

        # create triangular filters for each band
        bw = diff(band_edges)
        filterbank = zeros(nbands, length(sfreq))

        for k in 1:nbands
            # rising side of triangle
            @. filterbank[k, p[k]:(p[k + 1] - 1)] = (sfreq[p[k]:(p[k + 1] - 1)] - band_edges[k]) / bw[k]
            # falling side of triangle
            @. filterbank[k, p[k + 1]:(p[k + 2] - 1)] = (band_edges[k + 2] - sfreq[p[k + 1]:(p[k + 2] - 1)]) / bw[k + 1]
        end

        bw = (band_edges[3:end] - band_edges[1:(end - 2)])

        # normalization
        norm != :none && normalize!(filterbank, norm, bw)
    end
    
    MelFb(sr, MelFbSetup(nbands, scale, norm, freqrange, semitonerange), MelFbData(filterbank, filter_freq))
end

function Base.show(io::IO, mel_fb::MelFb)
    print(io, "MelFbank(")
    print(io, "sr=$(mel_fb.sr), ")
    print(io, "nbands=$(mel_fb.setup.nbands), ")
    print(io, "scale=:$(mel_fb.setup.scale), ")
    print(io, "norm=:$(mel_fb.setup.norm), ")
    print(io, "freqrange=$(mel_fb.setup.freqrange), ")
    print(io, "fbank=$(size(mel_fb.data.fbank)), ")
    print(io, "freq=$(length(mel_fb.data.freq)) frequencies)")
end

function Base.display(mel_fb::MelFb)
    f_length = size(mel_fb.data.fbank, 2)
    freqs = range(0, round(Int, mel_fb.sr / 2), length=f_length)

    p = plot(title = "Filter Bank Responses", xlabel = "Frequency (Hz)", ylabel = "Amplitude", legend = false)

    for i in eachrow(mel_fb.data.fbank)
        plot!(p, freqs, i[1:f_length], label = "", alpha = 0.6)
    end

    display(p)
end

function get_melfb(; source::Stft, kwargs...)
    _get_melfb(;
        sr = source.sr, 
        source = vec(sum(source.data.spec, dims=2)), 
        sfreq = source.data.freq, 
        nfft = source.setup.nfft, 
        kwargs...
    )
end