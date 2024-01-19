# using SpecialFunctions
# using Statistics
# using FFTW

"""
Continuous 1-D wavelet transform

CWTFILTERBANK methods:

wt                  - Continuous wavelet transform
freqz               - Wavelet frequency responses
timeSpectrum        - Time-averaged wavelet spectrum
scaleSpectrum       - Scale-averaged wavelet spectrum
wavelets            - Time-domain wavelets
scales              - Wavelet scales
waveletsupport      - Wavelet time support
qfactor             - Wavelet Q-factor
powerbw             - 3-dB bandwidths of wavelet bandpass filters
centerFrequencies   - Wavelet bandpass center frequencies
centerPeriods       - Wavelet bandpass center periods

CWTFILTERBANK properties:

SamplingFrequency   - Sampling frequency
SignalLength        - Signal length
Wavelet             - Analysis wavelet
FrequencyLimits    -  Frequency limits (tuple)
VoicesPerOctave     - Voices per octave
Boundary            - Reflect or treat data as periodic
da fare:
SamplingPeriod      - Sampling period
PeriodLimits        - Period limits
TimeBandwidth       - Time-bandwidth product
WaveletParameters   - Morse wavelet parameters
"""

if (!@isdefined(fbCell))
    struct fbCell
        sr::Int64
        length::Int64
        wavelet::Symbol
        frqLimits::Tuple{Int64,Int64}
        vpo::Int64 # voices per octave
        boundary::Symbol
    end
end

if (!@isdefined(fbParameters))
    struct fbParameters
        wavelet::Symbol
        waveletParameters::AbstractArray{Any}
        length::Int64
        sr::Int64
        samplingPeriod::AbstractArray{Any}
        vpo::Int64
        timeBandwidth::Int64
        freqLimits::AbstractArray{Any}
        periodLimits::AbstractArray{Any}
        boundary::Symbol
    end
end

if (!@isdefined(fbData))
    struct fbData{T<:Real,S<:Real}
        ga::Int64
        be::Int64
        cutOff::Int64
        normfreqflag::Bool
        signalPad::Int64
        waveletCF::Real
        omega::AbstractVector{T}
        frequencies::AbstractVector{T}
        scales::AbstractVector{T}
        psiDFT::AbstractMatrix{S}
        waveletCenterFrequencies::AbstractVector{T}
        nyquistBin::Int64
        sigvar::Real
        # npad
        # PhiDFT
        # PsiHalfPowerBandwidth
        # PsiHalfPowerFrequencies
        # PlotString 
        # CurrentClass
    end
end

function fzero(
    funFcn::T,
    x::Tuple{S,S}
) where {T<:Function,S<:Real} # Single-variable nonlinear zero finding

    # initialization
    tol = eps(Float64)

    a, b = x
    c = NaN
    d = NaN
    e = NaN
    fa = funFcn(a)
    fb = funFcn(b)
    fc = fb

    # main loop, exit from middle of the loop
    while ((fb != 0) && (a != b))
        # insure that b is the best result so far, a is the previous value of b, and c is on the opposite side of the zero from b.
        if ((fb > 0) == (fc > 0))
            c = a
            fc = fa
            d = b - a
            e = d
        end
        if (abs(fc) < abs(fb))
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end

        # convergence test and possible exit
        m = 0.5 * (c - b)
        toler = 2.0 * tol * max(abs(b), 1.0)
        if ((abs(m) <= toler) || (fb == 0.0))
            return b
        end

        # choose bisection or interpolation
        if ((abs(e) < toler) || (abs(fa) <= abs(fb)))
            # bisection
            d = m
            e = m
        else
            # interpolation
            s = fb / fa
            if (a == c)
                # linear interpolation
                p = 2.0 * m * s
                q = 1.0 - s
            else
                # inverse quadratic interpolation
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            end
            if (p > 0)
                q = -q
            else
                p = -p
            end
            # is interpolated point acceptable
            if ((2.0 * p < 3.0 * m * q - abs(toler * q)) && (p < abs(0.5 * e * q)))
                e = d
                d = p / q
            else
                d = m
                e = m
            end
        end # interpolation

        # next point
        a = b
        fa = fb
        if (abs(d) > toler)
            b = b + d
        elseif b > c
            b = b - toler
        else
            b = b + toler
        end
        fb = funFcn(b)
    end # main loop

    return b
end

function getFreqFromCutoffMorse(
    cutoff::T,
    cf::T,
    ga::Int64,
    be::Int64
) where {T<:Real}

    anorm = 2 * exp(be / ga * (1 + (log(ga) - log(be))))
    alpha = 2 * cutoff
    omax = ((750)^(1 / ga))

    psihat(om) = alpha - anorm * om^be * exp(-om^ga)

    if (psihat(cf) >= 0)
        if (psihat(omax) == psihat(cf))
            omegac = omax
        else
            omegac = cf
        end
    else
        omegac = fzero(psihat, (cf, omax))
    end
end # function getFreqFromCutoffMorse

function cwtfreqlimits(
    wavelet::Symbol,
    signalLength::Int64,
    ga::Int64,
    be::Int64,
    vpo::Int64,
    p::Int64,
    cutoff::Int64,
    fourierFactor::T,
    sigmaT::T,
    cf::T
) where {T<:Real}

    t = 1 #seconds
    fs = 1
    cutoff = cutoff / 100
    maxscale = signalLength / (sigmaT * p)

    if (wavelet == :morse)
        omegac = getFreqFromCutoffMorse(cutoff, cf, ga, be)
    elseif (wavelet == :bump)
        omegac = getFreqFromCutoffBump(cutoff, cf)
    elseif (wavelet == :amor)
        omegac = getFreqFromCutoffAmor(cutoff, cf)
    else
        error("Unknown wavelet ", wavelet, ".")
    end

    minscale = omegac / pi

    if (maxscale < minscale * 2^(1 / vpo))
        maxscale = minscale * 2^(1 / vpo)
    end

    minperiod = minscale * fourierFactor * t
    maxfreq = 1 / (minscale * fourierFactor) * fs

    maxperiod = maxscale * fourierFactor * t
    minfreq = 1 / (maxscale * fourierFactor) * fs

    if ((maxfreq > fs / 2) || (minperiod < 2 * t))
        maxfreq = fs / 2
        minperiod = 2 * t
    end

    return minfreq, maxperiod, maxscale, minscale, maxfreq, minperiod
end # function cwtfreqlimits

function wavCFandSD(
    wName::Symbol,
    ga::Int64,
    be::Int64
)
    cf = 0
    sigmaT = 0

    if (wName == :morse)
        cf = exp(1 / ga * (log(be) - log(ga)))

        # da morseproperties
        frac(a, b) = a / b
        morse_loga(a, b) = frac(b, a) .* (1 + log(a) - log(b))

        logsigo1 = frac(2, ga) .* log(frac(ga, 2 * be)) + loggamma(frac(2 * be + 1 + 2, ga)) - loggamma(frac(2 * be + 1, ga))
        logsigo2 = frac(2, ga) .* log(frac(ga, 2 * be)) + 2 .* loggamma(frac(2 * be + 2, ga)) - 2 .* loggamma(frac(2 * be + 1, ga))

        sigo = sqrt(exp(logsigo1) - exp(logsigo2))
        ra = 2 * morse_loga(ga, be) - 2 * morse_loga(ga, be - 1) + morse_loga(ga, 2 * (be - 1)) - morse_loga(ga, 2 * be)
        rb = 2 * morse_loga(ga, be) - 2 * morse_loga(ga, be - 1 + ga) + morse_loga(ga, 2 * (be - 1 + ga)) - morse_loga(ga, 2 * be)
        rc = 2 * morse_loga(ga, be) - 2 * morse_loga(ga, be - 1 + ga ./ 2) + morse_loga(ga, 2 * (be - 1 + ga ./ 2)) - morse_loga(ga, 2 * be)

        logsig2a = ra + frac(2, ga) .* log(frac(be, ga)) + 2 * log(be) + loggamma(frac(2 * (be - 1) + 1, ga)) - loggamma(frac(2 * be + 1, ga))
        logsig2b = rb + frac(2, ga) .* log(frac(be, ga)) + 2 * log(ga) + loggamma(frac(2 * (be - 1 + ga) + 1, ga)) - loggamma(frac(2 * be + 1, ga))
        logsig2c = rc + frac(2, ga) .* log(frac(be, ga)) + log(2) + log(be) + log(ga) + loggamma(frac(2 * (be - 1 + ga ./ 2) + 1, ga)) - loggamma(frac(2 * be + 1, ga))

        sig2a = exp(logsig2a)
        sig2b = exp(logsig2b)
        sig2c = exp(logsig2c)
        sigt = sqrt(sig2a + sig2b - sig2c)

        sigmaT = Real(sigt)
    elseif (wName == :amor)
        cf = 6
        sigmaT = sqrt(2)
    else
        cf = 5
        sigmaT = 5.847705
    end

    fourierFactor = (2 * pi) / cf

    return fourierFactor, sigmaT, cf
end # function wavCFandSD

function wfilterbank(
    fbcell::fbCell,
    scales::AbstractVector{T},
    omega::AbstractVector{T},
    ga::Int64,
    be::Int64) where {T<:Real} #Wavelet Filter Bank

    if (fbcell.wavelet == :morse)
        somega = scales * omega'
        absomega = abs.(somega)
        powscales = absomega .^ ga

        peakAF = exp(1 / ga * (log(be) - log(ga)))
        peakCF = peakAF / (2 * pi)

        factor = exp(-be * log(peakAF) + peakAF^ga)
        psidft = 2 * factor * exp.(be * log.(absomega) .- powscales) .* (somega .> 0)

        f = (peakAF ./ scales) / (2 * pi)
    else

    end

    return psidft, f
end # function wfilterbank

function cwtfilterbank(
    fbcell::fbCell,
    ga::Int64,
    be::Int64)

    # setup parameters
    timeBandwidth = ga * be
    p = 2 # Number of standard deviations

    fb_parameters = fbParameters(fbcell.wavelet, [], fbcell.length, sr, [], fbcell.vpo, timeBandwidth, [], [], fbcell.boundary)

    if (fbcell.wavelet == :morse)
        cutoff = 50
    else
        cutoff = 10
    end

    normfreqflag = false

    if (fbcell.boundary == :reflection)
        if fbcell.length <= 1e5
            signalPad = Int64(floor(fbcell.length / 2))
        else
            signalPad = Int64(ceil(log2(fbcell.length)))
        end
    else
        signalPad = 0
    end

    fourierFactor, sigmaT, cf = wavCFandSD(fbcell.wavelet, ga, be)

    # frequency grid
    N = fbcell.length + 2 * signalPad

    omega = [1:Int(floor(N / 2))...]
    omega = omega .* ((2 * pi) / N)
    omega = vcat(0.0, omega, -omega[Int(floor((N - 1) / 2)):-1:1])

    frequencies = fbcell.sr * omega ./ (2 * pi)

    minFreq, maxPeriod, maxScale, minScale, maxFreq, minPeriod = cwtfreqlimits(fbcell.wavelet, fbcell.length, ga, be, fbcell.vpo, p, cutoff, fourierFactor, sigmaT, cf)
    numoctaves = max(log2(maxScale / minScale), 1 / fbcell.vpo)
    a0 = 2^(1 / fbcell.vpo)
    scales = minScale * a0 .^ (0:numoctaves*fbcell.vpo)

    # compute filter bank
    psidft, f = wfilterbank(fbcell, scales, omega, ga, be)

    # nyquistBin
    nt = round(Int, size(psidft, 2))
    f = f .* fbcell.sr
    nyquistBin = (nt >>> 1) + 1

    return (; ga, be, cutoff, normfreqflag, signalPad, cf, omega, frequencies, scales, psidft, f, nyquistBin)

end # function cwtfilterbank

function createCoiIndices(
    n::Int64
)

    indices = vec(zeros(n, 1))
    if (isodd(length(x)))
        # odd length case
        M = Int64(ceil(n / 2))
        indices[1:M] .= [1:M;] #[] e ; genera una serie di numeri
        indices[M+1:N] .= [M-1:-1:1;]
    else
        # even length case
        indices[1:Int64(n / 2)] = [1:Int64(n / 2);]
        indices[Int64(n / 2)+1:n] .= [Int64(n / 2):-1:1;]
    end

    return indices
end # function createCoiIndices

function wt(
    x::Union{AbstractVector{T},AbstractArray{T}},
    fbcell::fbCell,
    fbdataT::NamedTuple
) where {T<:Real}

    n = fbcell.length
    dataclass = eltype(x)
    psihat = fbdataT.psidft
    # check whether input is real or complex
    isRealX = isreal(x)

    # x = hcat(x...) #trasforma da vettore a matrice ad una riga
    sigvar = var(x, corrected=false) # per avere lo stesso risultato di matlab

    xv = x
    if (fbdataT.signalPad > 0)
        xv = vcat(reverse(xv[1:fbdataT.signalPad]), xv, xv[end:-1:end-fbdataT.signalPad+1])
    end

    # fourier transform of input
    xposdft = fft(xv)
    xposdft = hcat(xposdft...)
    # obtain the CWT in the Fourier domain
    cfsposdft = xposdft .* fbdataT.psidft
    # invert to obtain wavelet coefficients
    cfspos = ifft(cfsposdft, 2)
    cfs = cfspos

    if (fbdataT.signalPad > 0)
        cfs = cfs[:, fbdataT.signalPad+1:fbdataT.signalPad+n, :]
    end

    f = fbdataT.f
    if (fbcell.wavelet == :morse)
        FourierFactor, sigmaPsi = wavCFandSD(fbcell.wavelet, fbdataT.ga, fbdataT.be)
    else

    end
    coiScalar = FourierFactor / sigmaPsi

    dt = 1 / fbcell.sr
    samples = createCoiIndices(n)
    coitmp = coiScalar * dt * samples
    coi = 1 ./ coitmp
    max_coi = max(fbdataT.f...)
    for i in eachindex(coi)
        if coi[i] > max_coi
            coi[i] = max_coi
        end
    end

    return sigvar, cfs, f, coi
end # function wt

function cwt(
    x::Union{AbstractVector{T},AbstractArray{T}},
    sr::Int64,
    wavelet::Symbol=:morse,
    sigLen::Int64=1024,
    vpo::Int64=10, # VoicesPerOctave
    boundary::Symbol=:reflection,
    frqLimitLow::Int64=0,
    frqLimitHi::Int64=Int(round(sr / 2)),
    ga::Int64=3,
    be::Int64=20
) where {T<:Real}

    signalLength = size(x, 1) # lavora solo con file mono
    fbcell = fbCell(sr, signalLength, wavelet, (frqLimitLow, frqLimitHi), vpo, boundary)
    fbdataT = cwtfilterbank(fbcell, ga, be)

    # sigvar, cfs, freq, coitmp = wt(x, fbcell, fbdataT)
    sigvar, cfs, freq, coitmp = wt(x, fbcell, fbdataT)

    # dt = 1/sr
    # t = [0:dt:signalLength*dt-dt;]

    # fbdata = fbData(fbdataT..., sigvar)

    return cfs, freq, coitmp #,scalcfs
end # function cwt

# # debug
# using PyCall
# librosa = pyimport("librosa")
# sr_src = 8000
# x, sr = librosa.load("/home/riccardopasini/Documents/Aclai/Julia_additional_files/test.wav", sr=sr_src, mono=true)

# depth = 5
# wname = "db2"
# spec, freqs, times = cwt(x, sr)