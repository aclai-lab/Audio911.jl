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
    peakCF = peakAF / (2 * pi)

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
    alpha = 2 * cutoff / 100

    psihat = x -> alpha - 2 * exp(-(x - cf) .^ 2 / 2)

    omax = ((2 * 750) .^ 0.5 + cf)

    if psihat(cf) > 0
        omegaC = omax
    else
        omegaC = find_zero(psihat, (cf, omax))
    end
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
function frequency_grid(
    sr::Int64,
    x_length::Int64,
    signal_pad::Int64
)
    n = x_length + 2 * signal_pad

    omega = [1:floor(Int, n / 2)...] .* ((2 * pi) / n)
    omega = vcat(0.0, omega, -omega[floor(Int, (n - 1) / 2):-1:1])
    frequencies = sr * omega ./ (2 * pi)

    return omega, frequencies
end

function cwt_scales(
    wavelet::Symbol,
    x_length::Int64,   
    gamma_beta::Tuple{Real, Real},
    vpo::Int64
)
    cutoff = wavelet == :morse ? 50 : 10
    ga, be = gamma_beta

    if wavelet == :morse
        center_freq, _ = morsepeakfreq(ga, be)
        _, _, _, sigmaT, _ = morseproperties(ga, be)

        omegaC = get_freq_cutoff_morse(cutoff, center_freq, ga, be)

    elseif wavelet == :morlet
        center_freq = 6
        sigmaT = sqrt(2)

        omegaC = get_freq_cutoff_morlet(cutoff, center_freq)

    elseif wavelet == :bump
        center_freq = 5
        # measured standard deviation of bump wavelet
        sigmaT = 5.847705

        omegaC = get_freq_cutoff_bump(cutoff, center_freq)
    else
        error("Wavelet $wavelet not supported.")
    end

    maxscale = x_length / (sigmaT * 2)
    minscale = omegaC / pi

    # if the max scale (min freq) is beyond the max freq, set it one step away
    if maxscale < minscale * 2^(1 / vpo)
        maxscale = minscale * 2^(1 / vpo)
    end

    numoctaves = log2(maxscale / minscale)
    a0 = 2^(1 / vpo)

    return minscale * a0 .^ (0:(vpo * numoctaves)), center_freq
end

# ---------------------------------------------------------------------------- #
#                       continuous wavelets filterbank                         #
# ---------------------------------------------------------------------------- #
function _get_cwt_fb(
    sr::Int64,
    x_length::Int64;
    wavelet::Symbol,
    morse_params::Tuple{Int64, Int64},
    vpo::Int64,
    boundary::Symbol
)
    gamma_beta = (morse_params[1], morse_params[2] / morse_params[1])
    ga, be = gamma_beta
    signal_pad = boundary == :reflection ? x_length <= 1e5 ? floor(Int, x_length / 2) : ceil(Int, log2(x_length)) : 0

    omega, frequencies = frequency_grid(sr, x_length, signal_pad)
    scales, center_freq = cwt_scales(wavelet, x_length, gamma_beta, vpo)

    somega = scales .* omega'

    if wavelet == :morse
        absomega = abs.(somega)
        if ga == 3
            powscales = absomega .* absomega .* absomega
        else
            powscales = absomega .^ ga
        end
        factor = exp(-be * log(center_freq) + center_freq^ga)
        fbank = 2 * factor * exp.(be .* log.(absomega) - powscales) .* (somega .> 0)

    elseif wavelet == :morlet
        fc = 6
        mul = 2
        squareterm = (somega .- fc) .* (somega .- fc)
        gaussexp = -squareterm ./ 2
        expnt = gaussexp .* (somega .> 0)
        fbank = mul * exp.(expnt) .* (somega .> 0)

    elseif wavelet == :bump
        fc = 5
        sigma = 0.6
        w = (somega .- fc) ./ sigma
        absw2 = w .* w
        expnt = -1 ./ (1 .- absw2)
        daughter = 2 * exp(1) * exp.(expnt) .* (abs.(w) .< 1 .- eps(1.0))
        daughter[isnan.(daughter)] .= 0
        fbank = daughter

    else
        error("Wavelet $wavelet not supported.")
    end

    wavelet_center_freqs = (center_freq ./ scales) / (2 * pi) .* sr

    return fbank, wavelet_center_freqs, signal_pad
end

# ---------------------------------------------------------------------------- #
#                                  callings                                    #
# ---------------------------------------------------------------------------- #
function get_cwt_fb!(
    rack::AudioRack,
    audio::Audio;
    wavelet::Symbol = :morse,
    morse_params::Tuple{Int64, Int64} = (3, 60),
    vpo::Int64 = 10,
    boundary::Symbol = :reflection
)
    fbank, wcf, signal_pad = _get_cwt_fb(rack.audio.sr, size(rack.audio.data, 1); wavelet, morse_params, vpo, boundary)

    rack.cwt_fb = CwtFbank(
        wavelet,
        morse_params,
        vpo,
        boundary, 
        signal_pad,
        fbank, 
        reverse(wcf)
    )
end

# TODO: calling for stft