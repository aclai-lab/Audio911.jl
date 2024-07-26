# ---------------------------------------------------------------------------- #
#                                cwt filterbank                                #
# ---------------------------------------------------------------------------- #

mutable struct CwtFbank
	# setup
    sr::Int64
	wavelet::Symbol # :morse, :morlet, :bump
	morse_params::Tuple{Int64, Int64}
	vpo::Int64
    freq_range::Tuple{Int64, Int64}
	# data
	fbank::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function CwtFbank(;
    sr,
	wavelet = :morse,
	morse_params = (3, 60),
	vpo = 10,
    freq_range = (20, round(Int, sr / 2)),
    fbank = nothing,
    freq = nothing,
)
    cwt_fb = CwtFbank(
        sr,
        wavelet,
        morse_params,
        vpo,
        freq_range,
        fbank,
        freq
    )
	return 	cwt_fb
end

# the wavelets are normalized so that the peak magnitudes for all passbands are approximately equal to 2.

# vpo = voices per octave

# morse wavelet parameters, specified as a two-element vector. 
# the first element is the symmetry parameter (gamma), which must be greater than or equal to 1. 
# the second element is the time-bandwidth product, which must be greater than or equal to gamma. 
# the ratio of the time-bandwidth product to gamma cannot exceed 40.
# when gamma is equal to 3, the Morse wavelet is perfectly symmetric in the frequency domain. 
# the skewness is equal to 0. Values of gamma greater than 3 result in positive skewness, 
# while values of gamma less than 3 result in negative skewness.

# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
function morsepeakfreq(ga::Real, be::Real)
    # peak frequency for 0-th order Morse wavelet is $(\frac{\beta}{\gamma})^{1/\gamma}$
    peakAF = exp(1 / ga * (log(be) - log(ga)))
    # obtain the peak frequency in cyclical frequency
    peakCF = peakAF / 2π

    return peakAF, peakCF
end

function morseproperties(ga::Real, be::Real)
    width = sqrt(ga * be)
    skew = (ga - 3) / width
    kurt = 3 - skew .^ 2 - (2 / width^2)
    
    lg, lb = log(ga), log(be)
    morse_loga = (be) -> be / ga * (1 + lg - log(be))
    morse_gb = 2 * morse_loga(be) - morse_loga(2 * be)

    lg_sig1 = 2 / ga * log(ga / (2 * be))
    lg_sig2 = 2 / ga * log(be / ga)
    lg_sig3 = loggamma((2 * be + 1) / ga)

    logsigo1 = lg_sig1 + loggamma((2 * be + 1 + 2) / ga) - lg_sig3
    logsigo2 = lg_sig1 + 2 * loggamma((2 * be + 2) / ga) - 2 * lg_sig3

    sigmaF = sqrt(exp(logsigo1) - exp(logsigo2))

    ra = morse_gb - 2 * morse_loga(be - 1)           + morse_loga(2 * (be - 1))
    rb = morse_gb - 2 * morse_loga(be - 1 + ga)      + morse_loga(2 * (be - 1 + ga))
    rc = morse_gb - 2 * morse_loga(be - 1 + ga ./ 2) + morse_loga(2 * (be - 1 + ga ./ 2))

    logsig2a = ra + lg_sig2 + 2 * lb + loggamma((2 * (be - 1) + 1) / ga) - lg_sig3
    logsig2b = rb + lg_sig2 + 2 * lg + loggamma((2 * (be - 1 + ga) + 1) / ga) - lg_sig3
    logsig2c = rc + lg_sig2 + log(2) + lb + lg + loggamma((2 * (be - 1 + ga ./ 2) + 1) / ga) - lg_sig3

    sigmaT = sqrt(exp(logsig2a) + exp(logsig2b) - exp(logsig2c))

    return width, skew, kurt, sigmaT, sigmaF
end

function get_freq_cutoff_morse(cutoff::Int64, cf::Real, ga::Real, be::Real)
    anorm = 2 * exp(be / ga * (1 + (log(ga) - log(be))))
    alpha = 2 * (cutoff / 100)

    psihat = x -> alpha - anorm * x .^ be * exp(-x .^ ga)

    omax = ((750) .^ (1 / ga))
    if psihat(cf) >= 0
        if psihat(omax) == psihat(cf)
            omegaC = omax
        else
            omegaC = cf
        end
    else
        omegaC = find_zero(psihat, (cf, omax))
    end
end

function get_freq_cutoff_morlet(cutoff::Int64, cf::Real)
    alpha = 0.02 * cutoff
    omax = √1500 + cf
    psihat(x) = alpha - 2 * exp(-(x - cf)^2 / 2)
    
    psihat(cf) > 0 ? omax : find_zero(psihat, (cf, omax))
end

function get_freq_cutoff_bump(cutoff::Int64, cf::Real)
    sigma = 0.6
    cutoff = cutoff / 100

    if cutoff < 10 * eps(0.0)
        omegaC = cf + sigma - 10 * eps(cf + sigma)
    else
        alpha = 2 * cutoff

        psihat = x -> 1 / (1 - x^2) + log(alpha) - log(2) - 1

        epsilon = find_zero(psihat, (0 + eps(0.0), 1 - eps(1.0)))
        omegaC = sigma * epsilon + cf
    end
end

# ---------------------------------------------------------------------------- #
#              construct the frequency grid for the wavelet DFT                #
# ---------------------------------------------------------------------------- #
function cwt_scales(
    wavelet::Symbol,
    ga::Real,
    be::Real,
    vpo::Int64,
    freq_range::Tuple{Int64, Int64},
    sr::Int64,
    x_length::Int64
)
    cutoff = wavelet == :morse ? 50 : 10

    if wavelet == :morse
        center_freq, _ = morsepeakfreq(ga, be)
        _, _, _, sigmaT, _ = morseproperties(ga, be)

        omegaC = get_freq_cutoff_morse(cutoff, center_freq, ga, be)

    elseif wavelet == :morlet
        center_freq = 6
        sigmaT = √2

        omegaC = get_freq_cutoff_morlet(cutoff, center_freq)

    elseif wavelet == :bump
        center_freq = 5
        # measured standard deviation of bump wavelet
        sigmaT = 5.847705

        omegaC = get_freq_cutoff_bump(cutoff, center_freq)
    else
        error("Wavelet $wavelet not supported.")
    end

    minfreq = sigmaT * center_freq * sr / (π * x_length)

    if freq_range[1] < minfreq
        freq_range = (minfreq, freq_range[2])
    end

    wrange = freq_range .* (2π / sr)
    wfreq = (center_freq / wrange[2], center_freq / wrange[1])
    a0 = 2^(1 / vpo)
    n_octaves = log2(wfreq[2] / wfreq[1])
    scales = wfreq[1] * a0 .^ (0:(vpo * n_octaves))

    return scales, center_freq
end

# ---------------------------------------------------------------------------- #
#                       continuous wavelets filterbank                         #
# ---------------------------------------------------------------------------- #
function _get_cwtfb(;
    x_length::Int64,
    cwt_fb::CwtFbank,
)
    ga, be = (cwt_fb.morse_params[1], cwt_fb.morse_params[2] / cwt_fb.morse_params[1])

    omega = range(0, step=(2π / x_length), length=floor(Int, x_length / 2) + 1)
    scales, center_freq = cwt_scales(cwt_fb.wavelet, ga, be, cwt_fb.vpo, cwt_fb.freq_range, cwt_fb.sr, x_length)

    somega = scales .* omega'

    if cwt_fb.wavelet == :morse
        absomega = abs.(somega)
        powscales = ga == 3 ? absomega .^ 3 : absomega .^ ga
        factor = exp(-be * log(center_freq) + center_freq^ga)
        cwt_fb.fbank = 2 * factor * exp.(be .* log.(absomega) - powscales) .* (somega .> 0)

    elseif cwt_fb.wavelet == :morlet
        fc, mul = 6, 2
        squareterm = (somega .- fc) .^ 2
        expnt = -squareterm ./ 2 .* (somega .> 0)
        cwt_fb.fbank = mul * exp.(expnt) .* (somega .> 0)

    elseif cwt_fb.wavelet == :bump
        fc, sigma = 5, 0.6
        w = (somega .- fc) ./ sigma
        absw2 = w .^ 2
        expnt = -1 ./ (1 .- absw2)
        daughter = 2 * exp(1) * exp.(expnt) .* (abs.(w) .< 1 .- eps(1.0))
        daughter[isnan.(daughter)] .= 0
        cwt_fb.fbank = daughter

    else
        error("Wavelet $cwt_fb.wavelet not supported.")
    end

    cwt_fb.freq = (center_freq ./ scales) / (2π) .* cwt_fb.sr
    cwt_fb.fbank = hcat(cwt_fb.fbank, zeros(size(cwt_fb.fbank[:,2:end])))

    return cwt_fb
end

function Base.show(io::IO, cwt_fb::CwtFbank)
    print(io, "CwtFbank(")
    print(io, "sr=$(cwt_fb.sr), ")
    print(io, "wavelet=:$(cwt_fb.wavelet), ")
    print(io, "morse_params=$(cwt_fb.morse_params), ")
    print(io, "vpo=$(cwt_fb.vpo), ")
    print(io, "freq_range=$(cwt_fb.freq_range), ")
    print(io, "full_proc=$(cwt_fb.full_proc), ")
    if isnothing(cwt_fb.fbank)
        print(io, "fbank=nothing, ")
    else
        print(io, "fbank=$(size(cwt_fb.fbank)), ")
    end
    if isnothing(cwt_fb.freq)
        print(io, "freq=nothing)")
    else
        print(io, "freq=$(length(cwt_fb.freq)) frequencies)")
    end
end

function Base.display(cwt_fb::CwtFbank)
    if isempty(cwt_fb.fbank)
        println("Filter bank is empty.")
        return
    end

    f_length = round(Int, size(cwt_fb.fbank, 2) / 2)
    freqs = range(0, round(Int, cwt_fb.sr / 2), length=f_length)

    p = plot(title = "Filter Bank Responses", xlabel = "Frequency (Hz)", ylabel = "Amplitude", legend = false)

    for i in eachrow(cwt_fb.fbank)
        plot!(p, freqs, i[1:f_length], label = "", alpha = 0.6)
    end

    display(p)
end

function get_cwtfb(;
	audio::Audio,
	kwargs...
)
    _get_cwtfb(; x_length=size(audio.data, 1), cwt_fb=CwtFbank(; sr=audio.sr, kwargs...))
end