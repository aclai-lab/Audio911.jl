# ---------------------------------------------------------------------------- #
#                                mel filterbank                                #
# ---------------------------------------------------------------------------- #

mutable struct MelFbank
	# setup
    sr::Int64
    nbands::Int64
    scale::Symbol # :mel_htk, :mel_slaney, :erb, :bark
    norm::Symbol # :bandwidth, :area, :none
    freq_range::Tuple{Int64, Int64}
	# data
	fbank::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function MelFbank(;
    sr,
	nbands = 26,
    scale = :mel_htk,
    norm = :bandwidth,
    freq_range = (0, round(Int, sr / 2)),
    fbank = nothing,
    freq = nothing,
)
    mel_fb = MelFbank(
        sr,
        nbands,
        scale,
        norm,
        freq_range,
        fbank,
        freq
    )
	return 	mel_fb
end

# ---------------------------------------------------------------------------- #
#                         scale convertions functions                          #
# ---------------------------------------------------------------------------- #
function hz2mel(hz::Tuple{Int64, Int64}, style::Symbol)
    style == :mel_htk && return 2595 .* log10.(1 .+ hz ./ 700)
    style == :mel_slaney && begin
        lin_step = 200 / 3
        return @. ifelse(hz < 1000, hz / lin_step,
            log(hz / 1000) / (log(6.4) / 27) + (1000 / lin_step))
    end
    error("Unknown style ($style).")
end

function mel2hz(mel_range::Tuple{Float64, Float64}, n_bands::Int64, style::Symbol)
    mel = LinRange(mel_range[1], mel_range[2], n_bands + 2)
    style == :mel_htk && return hz = 700 * (exp10.(mel / 2595) .- 1)
    style == :mel_slaney && begin
        lin_step = 200 / 3
        cp_mel = 1000 / lin_step
        return @. ifelse(
            mel < cp_mel, mel * lin_step, 1000 * exp(log(6.4) / 27 * (mel - cp_mel)))
    end
    error("Unknown style ($style).")
end

function hz2erb(hz::Tuple{Int64, Int64})
    log(10) * 1000 / (24.673 * 4.368) .* log10.(1 .+ 0.004368 .* hz)
end

function erb2hz(erb_range::Tuple{Float64, Float64}, n_bands::Int64)
    erb = LinRange(erb_range[1], erb_range[2], n_bands)
    (10 .^ (erb / (log(10) * 1000 / (24.673 * 4.368))) .- 1) ./ 0.004368
end

function hz2bark(hz::Tuple{Int64, Int64})
    bark = 26.81 .* hz ./ (1960 .+ hz) .- 0.53
    bark = map(x -> x < 2 ? 0.85 * x + 0.3 : x > 20.1 ? 1.22 * x - 4.422 : x, bark)
end

function bark2hz(bark_range::Tuple{Float64, Float64}, n_bands::Int64)
    bark = LinRange(bark_range[1], bark_range[2], n_bands + 2)
    bark = map(x -> x < 2 ? (x - 0.3) / 0.85 : x > 20.1 ? (x + 0.22 * 20.1) / 1.22 : x, bark)
    hz = 1960 .* (bark .+ 0.53) ./ (26.28 .- bark)
end

# ---------------------------------------------------------------------------- #
#                                 normalization                                #
# ---------------------------------------------------------------------------- #
function normalize!(
        filterbank::AbstractArray{Float64},
        norm::Symbol,
        bw::AbstractVector{Float64}
)
    norm_funcs = Dict(
        :area => () -> sum(filterbank, dims = 2),
        :bandwidth => () -> bw / 2,
        :none => () -> 1
    )

    weight_per_band = get(norm_funcs, norm, norm_funcs[:none])()
    filterbank ./= (weight_per_band .+ (weight_per_band .== 0))
end

# ---------------------------------------------------------------------------- #
#                                   gammatone                                  #
# ---------------------------------------------------------------------------- #
function compute_gammatone_coeffs(sr::Int64, bands::AbstractVector{Float64})
    t = 1 / sr
    erb = bands ./ 9.26449 .+ 24.7
    filt = 1.019 * 2π * erb

    a0 = t
    a2 = 0
    b0 = 1

    b1 = -2 * cos.(2 * bands * π * t) ./ exp.(filt * t)
    b2 = exp.(-2 * filt * t)

    a11 = -(2 * t * cos.(2 * bands * π * t) ./ exp.(filt * t) +
            2 * sqrt(3 + 2^(3 / 2)) * t * sin.(2 * bands * π * t) ./ exp.(filt * t)) / 2
    a12 = -(2 * t * cos.(2 * bands * π * t) ./ exp.(filt * t) -
            2 * sqrt(3 + 2^(3 / 2)) * t * sin.(2 * bands * π * t) ./ exp.(filt * t)) / 2
    a13 = -(2 * t * cos.(2 * bands * π * t) ./ exp.(filt * t) +
            2 * sqrt(3 - 2^(3 / 2)) * t * sin.(2 * bands * π * t) ./ exp.(filt * t)) / 2
    a14 = -(2 * t * cos.(2 * bands * π * t) ./ exp.(filt * t) -
            2 * sqrt(3 - 2^(3 / 2)) * t * sin.(2 * bands * π * t) ./ exp.(filt * t)) / 2

    gain = abs.(
        (
        (-2 * exp.(4 * 1im * bands * π * t) .* t .+
         2 * exp.(-(filt .* t) .+ 2 * 1im * bands * π * t) .* t .*
         (cos.(2 * bands * π * t) .- sqrt(3 - 2^(3 / 2)) .*
                                     sin.(2 * bands * π * t))) .*
        (-2 * exp.(4 * 1im * bands * π * t) .* t .+
         2 * exp.(-(filt .* t) .+ 2 * 1im * bands * π * t) .* t .*
         (cos.(2 * bands * π * t) .+ sqrt(3 - 2^(3 / 2)) .*
                                     sin.(2 * bands * π * t))) .*
        (-2 * exp.(4 * 1im * bands * π * t) .* t .+
         2 * exp.(-(filt .* t) .+ 2 * 1im * bands * π * t) .* t .*
         (cos.(2 * bands * π * t) .-
          sqrt(3 + 2^(3 / 2)) .* sin.(2 * bands * π * t))) .*
        (-2 * exp.(4 * 1im * bands * π * t) .* t .+
         2 * exp.(-(filt .* t) .+ 2 * 1im * bands * π * t) .* t .*
         (cos.(2 * bands * π * t) .+ sqrt(3 + 2^(3 / 2)) .* sin.(2 * bands * π * t)))
    ) ./
        (-2 ./ exp.(2 * filt * t) .- 2 * exp.(4 * 1im * bands * π * t) .+
         2 * (1 .+ exp.(4 * 1im * bands * π * t)) ./ exp.(filt * t)) .^ 4
    )

    allfilts = ones(length(bands))
    fcoefs = [a0 * allfilts, a11, a12, a13, a14, a2 * allfilts, b0 * allfilts, b1, b2, gain]

    coeffs = zeros(4, 6, length(bands))

    a0 = fcoefs[1]
    a11 = fcoefs[2]
    a12 = fcoefs[3]
    a13 = fcoefs[4]
    a14 = fcoefs[5]
    a2 = fcoefs[6]
    b0 = fcoefs[7]
    b1 = fcoefs[8]
    b2 = fcoefs[9]
    gain = fcoefs[10]

    for ind in 1:length(bands)
        coeffs[:, :, ind] = [a0[ind]/gain[ind] a11[ind]/gain[ind] a2[ind]/gain[ind] b0[ind] b1[ind] b2[ind];
                             a0[ind] a12[ind] a2[ind] b0[ind] b1[ind] b2[ind];
                             a0[ind] a13[ind] a2[ind] b0[ind] b1[ind] b2[ind];
                             a0[ind] a14[ind] a2[ind] b0[ind] b1[ind] b2[ind]]
    end

    return coeffs
end

# ---------------------------------------------------------------------------- #
#                           design filterbank matrix                           #
# ---------------------------------------------------------------------------- #
function _get_melfb(;
        stft_length::Int64,
        freq::StepRangeLen{Float64},
        mel_fb::MelFbank
)
    if mel_fb.scale == :erb
        erb_range = hz2erb(mel_fb.freq_range)
        filter_freq = erb2hz(erb_range, mel_fb.nbands)
        coeffs = compute_gammatone_coeffs(mel_fb.sr, filter_freq)

        iirfreqz = (b, a, n) -> fft([b; zeros(n - length(b))]) ./ fft([a; zeros(n - length(a))])
        sosfilt = (c, n) -> reduce((x, y) -> x .* y, map(row -> iirfreqz(row[1:3], row[4:6], n), eachrow(c)))
        apply_sosfilt = (i) -> abs.(sosfilt(coeffs[:, :, i], stft_length))

        filterbank = hcat(map(apply_sosfilt, 1:mel_fb.nbands)...)'
        # Derive Gammatone filter bandwidths as a function of center frequencies
        bw = 1.019 * 24.7 * (0.00437 * filter_freq .+ 1)

        # normalization
        (mel_fb.norm != :none) && normalize!(filterbank, mel_fb.norm, bw)

        rem(stft_length, 2) == 0 ? filterbank[:, 2:(stft_length ÷ 2)] .*= 2 : filterbank[:, 2:(stft_length ÷ 2 + 1)] .*= 2
        filterbank = filterbank[:, 1:(stft_length ÷ 2 + 1)]

    elseif mel_fb.scale == :mel_htk || mel_fb.scale == :mel_slaney || mel_fb.scale == :bark
        (mel_fb.scale == :mel_htk || mel_fb.scale == :mel_slaney) ? begin
            mel_range = hz2mel(mel_fb.freq_range, mel_fb.scale)
            band_edges = mel2hz(mel_range, mel_fb.nbands, mel_fb.scale)
        end : begin
            bark_range = hz2bark(mel_fb.freq_range)
            band_edges = bark2hz(bark_range, mel_fb.nbands)
        end

        filter_freq = band_edges[2:(end - 1)]
        mel_fb.nbands = length(filter_freq)

        p = [findfirst(freq .> edge) for edge in band_edges]
        isnothing(p[end]) ? p[end] = length(freq) : nothing

        # create triangular filters for each band
        bw = diff(band_edges)
        filterbank = zeros(mel_fb.nbands, length(freq))

        for k in 1:mel_fb.nbands
            # rising side of triangle
            filterbank[k, p[k]:(p[k + 1] - 1)] .= (freq[p[k]:(p[k + 1] - 1)] .- band_edges[k]) ./ bw[k]

            # falling side of triangle
            filterbank[k, p[k + 1]:(p[k + 2] - 1)] .= (band_edges[k + 2] .- freq[p[k + 1]:(p[k + 2] - 1)]) ./ bw[k + 1]
        end

        bw = (band_edges[3:end] - band_edges[1:(end - 2)])

        # normalization
        (mel_fb.norm != :none) && normalize!(filterbank, mel_fb.norm, bw)

    else
        error("Unknown filterbank frequency scale ($mel_fb.scale).")
    end

    mel_fb.fbank = filterbank
    mel_fb.freq = filter_freq

    return mel_fb
end

function Base.show(io::IO, mel_fb::MelFbank)
    print(io, "MelFbank(")
    print(io, "sr=$(mel_fb.sr), ")
    print(io, "nbands=$(mel_fb.nbands), ")
    print(io, "scale=:$(mel_fb.scale), ")
    print(io, "norm=:$(mel_fb.norm), ")
    print(io, "freq_range=$(mel_fb.freq_range), ")
    if isnothing(mel_fb.fbank)
        print(io, "fbank=nothing, ")
    else
        print(io, "fbank=$(size(mel_fb.fbank)), ")
    end
    if isnothing(mel_fb.freq)
        print(io, "freq=nothing)")
    else
        print(io, "freq=$(length(mel_fb.freq)) frequencies)")
    end
end

function Base.display(mel_fb::MelFbank)
    if isempty(mel_fb.fbank)
        println("Filter bank is empty.")
        return
    end

    f_length = size(mel_fb.fbank, 2)
    freqs = range(0, round(Int, mel_fb.sr / 2), length=f_length)

    p = plot(title = "Filter Bank Responses", xlabel = "Frequency (Hz)", ylabel = "Amplitude", legend = false)

    for i in eachrow(mel_fb.fbank)
        plot!(p, freqs, i[1:f_length], label = "", alpha = 0.6)
    end

    display(p)
end

function get_melfb(;
    stft::Stft,
    kwargs...
)
    _get_melfb(; stft_length=stft.stft_length, freq=stft.freq, mel_fb=MelFbank(; sr=stft.sr, kwargs...))
end
