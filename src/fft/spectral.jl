# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
"""
    SpectralSetup <: AbstractSetup

Configuration structure for basic spectral features.

# Fields
- `sr::Int64`: Sample rate in Hz
"""
struct SpectralSetup <: AbstractSetup
    sr::Int64
end

"""
    SpectralFluxSetup <: AbstractSetup

Configuration structure for spectral flux computation.

# Fields
- `sr::Int64`: Sample rate in Hz
- `p::Int64`: Norm order for flux calculation (default: 2 for Euclidean distance)
"""
struct SpectralFluxSetup <: AbstractSetup
    sr::Int64
    p::Int64
end

"""
    SpectralRolloffSetup <: AbstractSetup

Configuration structure for spectral rolloff computation.

# Fields
- `sr::Int64`: Sample rate in Hz
- `threshold::Real`: Energy threshold as fraction of total energy (default: 0.95)
"""
struct SpectralRolloffSetup <: AbstractSetup
    sr::Int64
    threshold::Real
end

# ---------------------------------------------------------------------------- #
#                                    utils                                     #
# ---------------------------------------------------------------------------- #
# Internal utility functions for spectral computations
sum_s = (s) -> sum(s, dims=2)
sum_sfreq = (s, sfreq) -> sum(s .* sfreq', dims=2)
arithmetic_mean = (s) -> sum_s(s) ./ size(s, 2)
centroid = (s, sfreq) -> vec(sum_sfreq(s, sfreq) ./ sum_s(s))
higher_moment_tmp = (s, sfreq) -> sfreq .- centroid(s, sfreq)'
spread = (s, sfreq) -> vec(sqrt.(sum((higher_moment_tmp(s, sfreq) .^ 2) .* s', dims=1) ./ sum_s(s)'))
higher_moment_denom = (s, sfreq) -> (spread(s, sfreq) .^ 3) .* sum_s(s)
higher_momement_num = (s, sfreq) -> (higher_moment_tmp(s, sfreq) .^ 3) .* s'

# ---------------------------------------------------------------------------- #
#                               spectral centroid                              #
# ---------------------------------------------------------------------------- #
"""
    SpectralCentroid{T} <: AbstractSpectral

Spectral centroid of an audio signal.

The spectral centroid indicates the center of mass of the spectrum. It is computed as 
the weighted mean of the frequencies present in the signal, weighted by their magnitudes.
Higher centroid values indicate brighter sounds with more high-frequency content.

# Fields
- `spec::Vector{<:AbstractFloat}`: Centroid values per frame in Hz
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralCentroid(spec::AbstractSpectrogram)

Compute spectral centroid from a spectrogram.

# Formula
``C = \\frac{\\sum_{k} f_k \\cdot S_k}{\\sum_{k} S_k}``

where `f_k` are frequency bins and `S_k` are magnitude values.
"""
struct SpectralCentroid{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralCentroid(
        spec::AbstractSpectrogram,
    )
        s, sfreq = get_data(spec), get_freq(spec)
        vspec = centroid(s, sfreq)
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                spectral crest                                #
# ---------------------------------------------------------------------------- #
"""
    SpectralCrest{T} <: AbstractSpectral

Spectral crest factor of an audio signal.

The spectral crest is the ratio of the maximum magnitude to the arithmetic mean of 
the magnitude spectrum. It measures the tonality of the signal: higher values indicate 
more tonal content with prominent peaks, while lower values suggest noise-like content.

# Fields
- `spec::Vector{<:AbstractFloat}`: Crest factor values per frame (dimensionless ratio)
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralCrest(spec::AbstractSpectrogram)

Compute spectral crest from a spectrogram.

# Formula
``\\text{Crest} = \\frac{\\max(S_k)}{\\frac{1}{N}\\sum_{k} S_k}``
"""
struct SpectralCrest{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralCrest(
        spec::AbstractSpectrogram,
    )
        s = get_data(spec)
        peak = maximum(s, dims=2)
        vspec = vec(peak ./ arithmetic_mean(s))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral decrease                              #
# ---------------------------------------------------------------------------- #
"""
    SpectralDecrease{T} <: AbstractSpectral

Spectral decrease of an audio signal.

The spectral decrease represents the amount of decrease in the spectral amplitude and
quantifies the steepness of the spectral slope. It gives more weight to lower frequencies.

# Fields
- `spec::Vector{<:AbstractFloat}`: Decrease values per frame
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralDecrease(spec::AbstractSpectrogram)

Compute spectral decrease from a spectrogram.

# Formula
``D = \\frac{1}{\\sum_{k=2}^{N} S_k} \\sum_{k=2}^{N} \\frac{S_k - S_1}{k-1}``

where `S_1` is the first bin magnitude.
"""
struct SpectralDecrease{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralDecrease(
        spec::AbstractSpectrogram,
    )
        s = get_data(spec)
        vspec = vec(sum((s[:, 2:end] .- s[:, 1]) ./ (1:size(s, 2)-1)', dims=2) ./ sum(s[:, 2:end], dims=2))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral entropy                               #
# ---------------------------------------------------------------------------- #
"""
    SpectralEntropy{T} <: AbstractSpectral

Spectral entropy of an audio signal.

Spectral entropy measures the disorder or randomness of the spectrum. Values range from 
0 to 1: low entropy indicates a tonal signal with distinct peaks, while high entropy 
indicates a noise-like signal with energy distributed across many frequencies.

# Fields
- `spec::Vector{<:AbstractFloat}`: Entropy values per frame (normalized, 0-1)
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralEntropy(spec::AbstractSpectrogram)

Compute spectral entropy from a spectrogram.

# Formula
``H = -\\frac{1}{\\log_2(N)} \\sum_{k} p_k \\log_2(p_k)``

where ``p_k = \\frac{S_k}{\\sum_j S_j}`` is the normalized spectrum.
"""
struct SpectralEntropy{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralEntropy(
        spec::AbstractSpectrogram,
    )
        s = get_data(spec)
        X = s ./ repeat(sum_s(s)', size(s, 2), 1)'
        t = replace!(-sum(X .* log2.(X), dims=2), NaN => 0)
        vspec = vec(t ./ log2(size(s, 2)))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral flatness                              #
# ---------------------------------------------------------------------------- #
"""
    SpectralFlatness{T} <: AbstractSpectral

Spectral flatness (Wiener entropy) of an audio signal.

Spectral flatness is the ratio of geometric mean to arithmetic mean of the spectrum.
Values range from 0 to 1: values near 0 indicate tonal signals, while values near 1 
indicate noise-like signals with flat spectra.

# Fields
- `spec::Vector{<:AbstractFloat}`: Flatness values per frame (0-1)
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralFlatness(spec::AbstractSpectrogram)

Compute spectral flatness from a spectrogram.

# Formula
``F = \\frac{\\exp(\\frac{1}{N}\\sum_{k} \\log S_k)}{\\frac{1}{N}\\sum_{k} S_k}``
"""
struct SpectralFlatness{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralFlatness(
        spec::AbstractSpectrogram,
    )
        s = get_data(spec)
        geometric_mean = exp.(sum(log.(s .+ eps(eltype(s))), dims=2) / size(s, 2))
        vspec = vec(geometric_mean ./ arithmetic_mean(s))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                 spectral flux                                #
# ---------------------------------------------------------------------------- #
"""
    SpectralFlux{T} <: AbstractSpectral

Spectral flux of an audio signal.

Spectral flux measures the rate of change of the power spectrum. It quantifies how 
quickly the spectral content changes over time. High flux values indicate rapid changes, 
often corresponding to note onsets or transients.

# Fields
- `spec::Vector{<:AbstractFloat}`: Flux values per frame
- `info::SpectralFluxSetup`: Configuration metadata including norm order

# Constructor
    SpectralFlux(spec::AbstractSpectrogram; p::Int64=2)

Compute spectral flux from a spectrogram.

# Arguments
- `spec::AbstractSpectrogram`: Input spectrogram
- `p::Int64=2`: Norm order (2 for Euclidean distance, 1 for Manhattan distance)

# Formula
``\\text{Flux}_t = ||S_t - S_{t-1}||_p``

where `S_t` is the spectrum at frame `t`.
"""
struct SpectralFlux{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralFluxSetup

    function SpectralFlux(
        spec::AbstractSpectrogram;
        p::Int64=2
    )
        s = get_data(spec)
        initial_condition = s[1, :]
        temp = diff(hcat(initial_condition, s'), dims=2)
        vspec = [LinearAlgebra.norm(temp[:, i], p) for i in axes(temp, 2)]
        info = SpectralFluxSetup(get_sr(spec), p)
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral kurtosis                              #
# ---------------------------------------------------------------------------- #
"""
    SpectralKurtosis{T} <: AbstractSpectral

Spectral kurtosis of an audio signal.

Spectral kurtosis measures the flatness or peakedness of the spectrum around its centroid.
It is the fourth normalized central moment. Positive kurtosis indicates a peaked spectrum,
while negative kurtosis indicates a flat spectrum. Normal distribution has kurtosis of 3.

# Fields
- `spec::Vector{<:AbstractFloat}`: Kurtosis values per frame
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralKurtosis(spec::AbstractSpectrogram)

Compute spectral kurtosis from a spectrogram.

# Formula
``K = \\frac{\\sum_{k} (f_k - C)^4 S_k}{(\\sigma^3 \\sum_{k} S_k) \\sigma}``

where `C` is the centroid and `σ` is the spread.
"""
struct SpectralKurtosis{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralKurtosis(
        spec::AbstractSpectrogram,
    )
        s, sfreq = get_data(spec), get_freq(spec)
        vspec = vec(sum(higher_momement_num(s, sfreq) .*
            higher_moment_tmp(s, sfreq), dims=1) ./
            (higher_moment_denom(s, sfreq) .* spread(s, sfreq))')
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral rolloff                               #
# ---------------------------------------------------------------------------- #
"""
    SpectralRolloff{T} <: AbstractSpectral

Spectral rolloff frequency of an audio signal.

The spectral rolloff is the frequency below which a specified percentage (default 95%) 
of the total spectral energy is contained. It provides a measure of spectral shape and 
is useful for distinguishing voiced from unvoiced speech and music from noise.

# Fields
- `spec::Vector{<:AbstractFloat}`: Rolloff frequencies per frame in Hz
- `info::SpectralRolloffSetup`: Configuration metadata including threshold

# Constructor
    SpectralRolloff(spec::AbstractSpectrogram; threshold::Real=0.95)

Compute spectral rolloff from a spectrogram.

# Arguments
- `spec::AbstractSpectrogram`: Input spectrogram
- `threshold::Real=0.95`: Energy threshold (0.85-0.95 typical, default 0.95 for 95%)
"""
struct SpectralRolloff{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralRolloffSetup

    function SpectralRolloff(
        spec::AbstractSpectrogram;
        threshold::Real=0.95
    )
        s, sfreq = get_data(spec), get_freq(spec)
        c = cumsum(s', dims=1)
        d = c[end, :] * threshold
        idx = [findfirst(real(c[:, i]) .>= real(d[i])) for i in eachindex(d)]
        vspec = sfreq[idx]
        info = SpectralRolloffSetup(get_sr(spec), threshold)
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                              spectral skewness                               #
# ---------------------------------------------------------------------------- #
"""
    SpectralSkewness{T} <: AbstractSpectral

Spectral skewness of an audio signal.

Spectral skewness measures the asymmetry of the spectrum around its centroid. It is 
the third normalized central moment. Positive skewness indicates energy concentrated 
below the centroid, while negative skewness indicates energy above the centroid.

# Fields
- `spec::Vector{<:AbstractFloat}`: Skewness values per frame
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralSkewness(spec::AbstractSpectrogram)

Compute spectral skewness from a spectrogram.

# Formula
``S = \\frac{\\sum_{k} (f_k - C)^3 S_k}{\\sigma^3 \\sum_{k} S_k}``

where `C` is the centroid and `σ` is the spread.
"""
struct SpectralSkewness{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralSkewness(
        spec::AbstractSpectrogram,
    )
        s, sfreq = get_data(spec), get_freq(spec)
        vspec = vec(sum(higher_momement_num(s, sfreq), dims=1) ./ higher_moment_denom(s, sfreq)')
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                               spectral slope                                 #
# ---------------------------------------------------------------------------- #
"""
    SpectralSlope{T} <: AbstractSpectral

Spectral slope of an audio signal.

The spectral slope is a linear approximation of the spectrum shape, computed by least 
squares regression. It represents the rate of change of the spectrum amplitude with 
frequency. Negative slopes indicate more energy in low frequencies.

# Fields
- `spec::Vector{<:AbstractFloat}`: Slope values per frame
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralSlope(spec::AbstractSpectrogram)

Compute spectral slope from a spectrogram.

# Formula
Slope is computed via least squares fit: ``m = \\frac{\\sum (S_k - \\bar{S})(f_k - \\bar{f})}{\\sum (f_k - \\bar{f})^2}``
"""
struct SpectralSlope{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralSlope(
        spec::AbstractSpectrogram,
    )
        s, sfreq = get_data(spec), get_freq(spec)
        f_minus_mu_f = sfreq .- sum(sfreq, dims=1) ./ size(s, 2)
        X_minus_mu_X = s .- arithmetic_mean(s)
        vspec = vec(eltype(s).(
            sum(X_minus_mu_X' .* f_minus_mu_f, dims=1) ./ sum(f_minus_mu_f .^ 2)))
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                spectral spread                               #
# ---------------------------------------------------------------------------- #
"""
    SpectralSpread{T} <: AbstractSpectral

Spectral spread (bandwidth) of an audio signal.

Spectral spread measures the concentration of the spectrum around its centroid. It is 
the second central moment (standard deviation) of the spectrum. Higher values indicate 
wider bandwidth, while lower values indicate concentrated spectral energy.

# Fields
- `spec::Vector{<:AbstractFloat}`: Spread values per frame in Hz
- `info::SpectralSetup`: Configuration metadata

# Constructor
    SpectralSpread(spec::AbstractSpectrogram)

Compute spectral spread from a spectrogram.

# Formula
``\\sigma = \\sqrt{\\frac{\\sum_{k} (f_k - C)^2 S_k}{\\sum_{k} S_k}}``

where `C` is the spectral centroid.
"""
struct SpectralSpread{T} <: AbstractSpectral
    spec::Vector{<:AbstractFloat}
    info::SpectralSetup

    function SpectralSpread(
        spec::AbstractSpectrogram,
    )
        s, sfreq = get_data(spec), get_freq(spec)
        vspec = spread(s, sfreq)
        info = SpectralSetup(get_sr(spec))
        new{typeof(spec)}(vspec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    get_data(s::AbstractSpectral) -> Vector

Get the spectral feature data vector.

Returns the computed spectral feature values for each frame of the input signal.
"""
@inline get_data(s::AbstractSpectral)  = s.spec

"""
    get_setup(s::AbstractSpectral) -> Union{SpectralSetup, SpectralFluxSetup, SpectralRolloffSetup}

Get the configuration metadata for the spectral feature.

Returns the setup structure containing sample rate and any feature-specific parameters.
"""
@inline get_setup(s::AbstractSpectral) = s.info
