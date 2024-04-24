# using SpecialFunctions
# using Statistics, Roots
# using FFTW
# using Parameters
# using Plots

"""
Continuous 1-D wavelet transform

TODO documentation
"""

###################################################################################
#                                data structures                                  #
###################################################################################

@with_kw mutable struct FbSetup
    sr::Int64
    length::Int64
    wavelet::Symbol = :morse # :morse, :amor, :bump
    gamma::Int64 = 3
    beta::Int64 = 20
    time_bandwidth::Int64 = 3 * 20
    cutoff::Int64 = 50
    vpo::Int64 = 10 # voices per octave
    boundary::Symbol = :reflection # :reflection, :periodic
    signal_pad::Int64 = 0
    frequency_range::Tuple{Int64, Int64} = (0, 0)
    center_freq::Float64 = 0.0
    omega::Vector{Float64} = []
    frequencies::Vector{Float64} = []
    scales::Vector{Float64} = []
    period_range::Tuple{Int64, Int64} = (0, 0) # TODO cwtfilterbank.m line: 1093
    sampling_range::AbstractArray{Float64} = [] # TODO if necessary
    psidft::AbstractArray{Float64} = []
    wavelet_center_freqs::AbstractArray{Float64} = []
end

function morsepeakfreq(ga::Int64, be::Int64)
    # peak frequency for 0-th order Morse wavelet is $(\frac{\beta}{\gamma})^{1/\gamma}$
    peakAF = exp(1 / ga * (log(be) - log(ga)))
    # obtain the peak frequency in cyclical frequency
    peakCF = peakAF / (2 * pi)

    return peakAF, peakCF
end

function morseproperties(ga::Int64, be::Int64)
    morse_loga = (ga, be) -> be / ga * (1 + log(ga) - log(be))

    width = sqrt(ga * be)
    skew = (ga - 3) / width
    kurt = 3 - skew .^ 2 - (2 / width^2)

    logsigo1 = 2 / ga * log(ga / (2 * be)) + loggamma((2 * be + 1 + 2) / ga) -
               loggamma((2 * be + 1) / ga)
    logsigo2 = 2 / ga * log(ga / (2 * be)) + 2 * loggamma((2 * be + 2) / ga) -
               2 * loggamma((2 * be + 1) / ga)

    sigmaF = sqrt(exp(logsigo1) - exp(logsigo2))
    ra = 2 * morse_loga(ga, be) - 2 * morse_loga(ga, be - 1) +
         morse_loga(ga, 2 * (be - 1)) - morse_loga(ga, 2 * be)
    rb = 2 * morse_loga(ga, be) - 2 * morse_loga(ga, be - 1 + ga) +
         morse_loga(ga, 2 * (be - 1 + ga)) - morse_loga(ga, 2 * be)
    rc = 2 * morse_loga(ga, be) - 2 * morse_loga(ga, be - 1 + ga ./ 2) +
         morse_loga(ga, 2 * (be - 1 + ga ./ 2)) - morse_loga(ga, 2 * be)

    logsig2a = ra + 2 / ga * log(be / ga) + 2 * log(be) +
               loggamma((2 * (be - 1) + 1) / ga) - loggamma((2 * be + 1) / ga)
    logsig2b = rb + 2 / ga * log(be / ga) + 2 * log(ga) +
               loggamma((2 * (be - 1 + ga) + 1) / ga) - loggamma((2 * be + 1) / ga)
    logsig2c = rc + 2 / ga * log(be / ga) + log(2) + log(be) + log(ga) +
               loggamma((2 * (be - 1 + ga ./ 2) + 1) / ga) - loggamma((2 * be + 1) / ga)

    sig2a = exp(logsig2a)
    sig2b = exp(logsig2b)
    sig2c = exp(logsig2c)
    sigmaT = sqrt(sig2a + sig2b - sig2c)

    return width, skew, kurt, sigmaT, sigmaF
end

function get_freq_cutoff_morse(cutoff::Int64, cf::Float64, ga::Int64, be::Int64)
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

function get_freq_cutoff_amor(cutoff::Int64, cf::Float64)
    alpha = 2 * cutoff

    psihat = x -> alpha - 2 * exp(-(x - cf) .^ 2 / 2)

    omax = ((2 * 750) .^ 0.5 + cf)

    if psihat(cf) > 0
        omegaC = omax
    else
        omegaC = find_zero(psihat, (cf, omax))
    end
end

function get_freq_cutoff_bump(cutoff::Int64, cf::Float64)
    sigma = 0.6

    if cutoff < 10 * eps(0.0)
        omegaC = cf + sigma - 10 * eps(cf + sigma)
    else
        alpha = 2 * cutoff

        psihat = x -> 1 / (1 - x^2) + log(alpha) - log(2) - 1

        epsilon = find_zero(psihat, (0 + eps(0.0), 1 - eps(1.0)))
        omegaC = sigma * epsilon + cf
    end
end

function freq2scales(
        sr::Int64, frequency_range::Tuple{Int64, Int64}, vpo::Int64, center_freq::Float64)
    # convert frequencies in Hz to radians/sample
    wrange = frequency_range .* (1 / sr * 2 .* pi)
    a0 = 2^(1 / vpo)
    s0 = center_freq / wrange[2]
    smax = center_freq / wrange[1]
    numoctaves = log2(smax / s0)
    scales = s0 * a0 .^ (0:(vpo * numoctaves))
end

function cwtfilterbank!(fb_setup::FbSetup)
    ###########################################################################
    #                             setup parameters                            #
    ###########################################################################
    fb_setup.cutoff = fb_setup.wavelet == :morse ? 50 : 10

    fb_setup.signal_pad = fb_setup.boundary == :reflection ?
                          fb_setup.length <= 1e5 ?
                          floor(Int, fb_setup.length / 2) :
                          ceil(Int, log2(fb_setup.length)) : 0

    if fb_setup.wavelet == :morse
        fb_setup.center_freq, _ = morsepeakfreq(fb_setup.gamma, fb_setup.beta)
        _, _, _, sigmaT, _ = morseproperties(fb_setup.gamma, fb_setup.beta)

        omegaC = get_freq_cutoff_morse(
            fb_setup.cutoff, fb_setup.center_freq, fb_setup.gamma, fb_setup.beta)

    elseif fb_setup.wavelet == :amor
        fb_setup.center_freq = 6
        sigmaT = sqrt(2)

        omegaC = get_freq_cutoff_amor(fb_setup.cutoff, fb_setup.center_freq)

    elseif fb_setup.wavelet == :bump
        fb_setup.center_freq = 5
        # measured standard deviation of bump wavelet
        sigmaT = 5.847705

        omegaC = get_freq_cutoff_bump(fb_setup.cutoff, fb_setup.center_freq)
    else
        omegaC = pi
    end

    maxscale = fb_setup.length / sigmaT / 2
    minscale = omegaC / pi
    # if the max scale (min freq) is beyond the max freq, set it one step away
    if maxscale < minscale * 2^(1 / fb_setup.vpo)
        maxscale = minscale * 2^(1 / fb_setup.vpo)
    end

    fourier_factor = (2 * pi) / fb_setup.center_freq
    t = 1 / fb_setup.sr

    minperiod = minscale * fourier_factor * t
    maxfreq = 1 / (minscale * fourier_factor) * fb_setup.sr

    maxperiod = maxscale * fourier_factor * t
    minfreq = 1 / (maxscale * fourier_factor) * fb_setup.sr

    # guard against edge case
    if maxfreq > fb_setup.sr / 2 || minperiod < 2 * t
        maxfreq = fb_setup.sr / 2
        minperiod = 2 * t
    end

    if fb_setup.frequency_range[1] < minfreq
        fb_setup.frequency_range[1] = minfreq
    end

    ###########################################################################
    #            construct the frequency grid for the wavelet DFT             #
    ###########################################################################
    n = fb_setup.length + 2 * fb_setup.signal_pad

    omega = [1:floor(Int, n / 2)...] .* ((2 * pi) / n)
    fb_setup.omega = vcat(0.0, omega, -omega[floor(Int, (n - 1) / 2):-1:1])
    fb_setup.frequencies = fb_setup.sr * fb_setup.omega ./ (2 * pi)

    fb_setup.scales = freq2scales(
        fb_setup.sr, fb_setup.frequency_range, fb_setup.vpo, fb_setup.center_freq)

    somega = fb_setup.scales * fb_setup.omega'

    if fb_setup.wavelet == :morse
        absomega = abs.(somega)
        if fb_setup.gamma == 3
            powscales = absomega .* absomega .* absomega
        else
            powscales = absomega .^ fb_setup.gamma
        end
        factor = exp(-fb_setup.beta * log(fb_setup.center_freq) +
                     fb_setup.center_freq^fb_setup.gamma)
        fb_setup.psidft = 2 * factor * exp.(fb_setup.beta .* log.(absomega) - powscales) .*
                          (somega .> 0)

    elseif fb_setup.wavelet == :amor
        fc = 6
        mul = 2
        squareterm = (somega .- fc) .* (somega .- fc)
        gaussexp = -squareterm ./ 2
        expnt = gaussexp .* (somega .> 0)
        fb_setup.psidft = mul * exp.(expnt) .* (somega .> 0)

    else
        fc = 5
        sigma = 0.6
        w = (somega .- fc) ./ sigma
        absw2 = w .* w
        expnt = -1 ./ (1 .- absw2)
        daughter = 2 * exp(1) * exp.(expnt) .* (abs.(w) .< 1 .- eps(1.0))
        daughter[isnan.(daughter)] .= 0
        fb_setup.psidft = daughter
    end

    f = (fb_setup.center_freq ./ fb_setup.scales) / (2 * pi)

    fb_setup.wavelet_center_freqs = f .* fb_setup.sr
    # nyquist_bin = (size(psidft, 2) >>> 1) + 1

end # function cwtfilterbank

function wt(
        x::AbstractVector{Float64},
        fb_setup::FbSetup
)
    if (fb_setup.signal_pad > 0)
        x = vcat(reverse(x[1:fb_setup.signal_pad]), x,
            x[end:-1:(end - fb_setup.signal_pad + 1)])
    end

    # fourier transform of input
    xposdft = fft(x)
    # obtain the CWT in the Fourier domain
    cfsposdft = xposdft' .* fb_setup.psidft
    # invert to obtain wavelet coefficients
    cfs = ifft(cfsposdft, 2)
    # cfs = cfspos

    if (fb_setup.signal_pad > 0)
        cfs[:, fb_setup.signal_pad + 1:fb_setup.signal_pad + fb_setup.length]
    end
end # function wt

function cwt_windowing(
    cwt_spectrum::Matrix{Float64},
    window_length::Int64
)
    c_length = size(cwt_spectrum, 1)
    n_feats = size(cwt_spectrum, 2)
    n_hops = floor(Int, c_length / window_length)

    y = zeros(Float64, n_hops, n_feats)

    for i = 1:n_feats
        for j = 1:n_hops

                # y[j,i] = maximum(cwt_spectrum[(j-1)*window_length+1:j*window_length, i])
                y[j,i] = mean(cwt_spectrum[(j-1)*window_length+1:j*window_length, i])

        end
    end

    return y
end # function buffer

function cwt(
        x::AbstractVector{Float64},
        sr::Int64;
        wavelet::Symbol = :morse,
        ga::Int64 = 3,
        be::Int64 = 20,
        frequency_range::Tuple{Int64, Int64} = (0, round(Int, sr / 2)),
        vpo::Int64 = 10, # VoicesPerOctave
        boundary::Symbol = :reflection
)
    # ga, be = gamma, beta: symmetric parameters for morse wavelet

    fb_setup = FbSetup(
        sr = sr,
        length = size(x, 1),
        wavelet = wavelet,
        gamma = ga,
        beta = be,
        time_bandwidth = ga * be,
        vpo = vpo,
        boundary = boundary,
        frequency_range = frequency_range
    )

    cwtfilterbank!(fb_setup)

    return wt(x, fb_setup), fb_setup.wavelet_center_freqs
end # function cwt

function cwt(x::AbstractVector{<:AbstractFloat}, sr::Int64; kwargs...)
    cwt(Float64.(x), sr; kwargs...)
end