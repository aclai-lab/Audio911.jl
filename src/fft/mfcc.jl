# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
"""
    MfccSetup <: AbstractSetup

Configuration structure for MFCC (and GtCC) computation.

# Fields
- `sr::Int64`: Sample rate of the source audio signal, in Hz.
- `ncoeffs::Int64`: Number of cepstral coefficients retained.
- `rect::Base.Callable`: Rectification function applied before DCT (`mlog` or `cubic_root`).
- `dither::Bool`: Whether dithering was applied to avoid `log(0)`.
"""
struct MfccSetup <: AbstractSetup
    sr      :: Int64
    ncoeffs :: Int64
	rect    :: Base.Callable
	dither  :: Bool
end

# ---------------------------------------------------------------------------- #
#                                     mfcc                                     #
# ---------------------------------------------------------------------------- #
"""
    Mfcc{T} <: AbstractCepstrum

Mel-Frequency Cepstral Coefficients (MFCCs) of an audio signal.

MFCCs are derived by applying a Discrete Cosine Transform (DCT) to the
log-compressed mel-frequency energies. They provide a compact representation
of the spectral envelope widely used in speech recognition, speaker
identification, and music analysis.

# Fields
- `spec::Matrix{T}`: Cepstral coefficient matrix `(ncoeffs × nframes)`.
- `info::MfccSetup`: Configuration metadata.

# See also
[`Gtcc`](@ref), [`MelSpec`](@ref), [`Delta`](@ref)
"""
struct Mfcc{T} <: AbstractCepstrum
    spec :: Matrix{T} 
    info :: MfccSetup

    function Mfcc(
        spec :: Matrix{T},
        info :: MfccSetup       
    ) where T
        new{T}(spec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                     gtcc                                     #
# ---------------------------------------------------------------------------- #
"""
    Gtcc{T} <: AbstractCepstrum

Gammatone Cepstral Coefficients (GTCCs) of an audio signal.

GTCCs are the cepstral counterpart of [`ErbSpec`](@ref): a DCT is applied to the
compressed ERB-scale gammatone spectrogram instead of a mel spectrogram. Because
the underlying gammatone filterbank models the auditory periphery more faithfully
than triangular mel filters, GTCCs can offer improved robustness in noisy
conditions and in biologically-inspired audio analysis tasks.

The computation pipeline mirrors that of [`Mfcc`](@ref) exactly; the only
difference is that the input must be an [`ErbSpec`](@ref).

# Fields
- `spec::Matrix{T}`: Cepstral coefficient matrix `(ncoeffs × nframes)`.
- `info::MfccSetup`: Configuration metadata (shared with [`Mfcc`](@ref)).

# Constructors

    Gtcc(erb::ErbSpec; ncoeffs, rect, dither)

Convenience constructor. Computes GTCCs directly from an [`ErbSpec`](@ref).

# See also
[`Mfcc`](@ref), [`ErbSpec`](@ref), [`Delta`](@ref)
"""
struct Gtcc{T} <: AbstractCepstrum
    spec :: Matrix{T}
    info :: MfccSetup

    function Gtcc(
        spec :: Matrix{T},
        info :: MfccSetup
    ) where T
        new{T}(spec, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    Base.eltype(::Mfcc{T}) -> Type
    Base.eltype(::Gtcc{T}) -> Type

Return the element type of the cepstrum spectrogram data.
"""
Base.eltype(::Mfcc{T}) where {T} = T
Base.eltype(::Gtcc{T}) where {T} = T

"""
    get_data(m::AbstractCepstrum) -> Matrix

Get the cepstrum spectrogram data matrix, transposed to (nframes × nbands).
"""
@inline get_data(m::AbstractCepstrum)  = m.spec'

"""
    get_setup(m::AbstractCepstrum) -> MelSpecSetup

Get the configuration metadata for the cepstrum spectrogram.
"""
@inline get_setup(m::AbstractCepstrum) = m.info

"""
    get_sr(m::AbstractCepstrum) -> Int64

Return the sample rate (in Hz) associated with the cepstrum spectrogram.
"""
@inline get_sr(m::AbstractCepstrum) = m.info.sr

# ---------------------------------------------------------------------------- #
#                                    utils                                     #
# ---------------------------------------------------------------------------- #
"""
    mlog(x, y) -> Matrix

Apply base-10 logarithmic scaling to spectral data for MFCC computation.

Performs element-wise logarithmic transformation: `x * log10.(y)`, where `x` is 
typically the DCT matrix and `y` is the mel-frequency spectrogram. This compresses 
the dynamic range and better approximates human auditory perception.
"""
mlog(x, y) =  x * log10.(y)

"""
    cubic_root -> Function

Apply cube root scaling to spectral data for MFCC computation.

Performs element-wise cube root transformation: `x * y .^ (1/3)`, where `x` is 
typically the DCT matrix and `y` is the mel-frequency spectrogram. Provides 
compression similar to logarithmic scaling but with different mathematical properties, 
sometimes preferred for robustness.
"""
cubic_root(x, y) = x * y .^ (1 / 3)


function _create_dct_matrix(mel_coeffs::Int64)
    matrix = zeros(Float64, mel_coeffs, mel_coeffs)
    s0  = sqrt(1 / mel_coeffs)
    s1  = sqrt(2 / mel_coeffs)
    pic = 2 * pi / (2 * mel_coeffs)

    matrix[1, :] .= s0
    for k in 1:mel_coeffs, n in 2:mel_coeffs
        matrix[n, k] = s1 * cos(pic * (n - 1) * (k - 0.5))
    end

    return matrix
end

# ---------------------------------------------------------------------------- #
#                                  _cepstrum                                   #
# ---------------------------------------------------------------------------- #
"""
    _cepstrum(spec::AbstractSpectrogram; ncoeffs, rect, dither) -> (Matrix, MfccSetup)

Internal function that computes cepstral coefficients from any spectrogram.

Applies dithering, builds the DCT matrix, applies the rectification function,
and returns the truncated coefficient matrix together with the setup struct.
Called by both [`Mfcc`](@ref) and [`Gtcc`](@ref).

# Arguments
- `spec::AbstractSpectrogram`: Input spectrogram (e.g. [`MelSpec`](@ref) or [`ErbSpec`](@ref)).

# Keywords
- `ncoeffs::Int64`: Number of cepstral coefficients to retain (default: `nbands ÷ 2`).
- `rect::Base.Callable=mlog`: Rectification function (`mlog` or `cubic_root`).
- `dither::Bool=false`: If `true`, values below `1e-8` are clamped to `1e-8`.
  If `false`, exact zeros are replaced with `floatmin(Float64)`.

# Returns
- `coeffs::Matrix`: Cepstral coefficient matrix `(ncoeffs × nframes)`.
- `info::MfccSetup`: Filled configuration struct.
"""
function _cepstrum(
    spec    :: AbstractSpectrogram;
    ncoeffs :: Int64=get_nbands(spec) ÷ 2,
    rect    :: Base.Callable=mlog,
    dither  :: Bool=false,
)
    s      = get_data(spec)       # (nframes × nbands)
    sr     = get_sr(spec)
    nbands = get_nbands(spec)

    dither ? s[s .< 1e-8] .= 1e-8 : s[s .== 0] .= floatmin(Float64)

    dct_matrix = _create_dct_matrix(nbands)
    coeffs     = rect(dct_matrix, s')   # (nbands × nframes)

    info = MfccSetup(sr, ncoeffs, rect, dither)

    return coeffs[1:ncoeffs, :], info
end

# ---------------------------------------------------------------------------- #
#                                   get mfcc                                   #
# ---------------------------------------------------------------------------- #
"""
    Mfcc(mel::MelSpec; ncoeffs, rect, dither) -> Mfcc

Compute Mel-Frequency Cepstral Coefficients (MFCCs) from a mel spectrogram.

# Arguments
- `mel::MelSpec`: Input mel spectrogram (must use `htk` or `slaney` scale).

# Keywords
- `ncoeffs::Int64`: Number of cepstral coefficients to retain
  (default: `nbands ÷ 2`).
- `rect::Base.Callable=mlog`: Rectification applied before DCT.
  - `mlog`: logarithmic compression (MATLAB default, HTK-compatible).
  - `cubic_root`: cube-root compression.
- `dither::Bool=false`: Clamp near-zero values to `1e-8` before rectification
  to avoid `log(0)`. When `false`, exact zeros are replaced with
  `floatmin(Float64)`.

# Examples
```julia
frames   = Frames(audio; winsize=512, winstep=256)
stft     = Stft(frames; spectrum=power)
mel      = MelSpec(stft; nbands=40, scale=htk)

# 13 MFCCs with log compression
mfcc = Mfcc(mel; ncoeffs=13)

# 20 MFCCs with cube-root compression and dithering
mfcc = Mfcc(mel; ncoeffs=20, rect=cubic_root, dither=true)
```

# See also
[`Gtcc`](@ref), [`MelSpec`](@ref), [`Delta`](@ref), [`mlog`](@ref)
"""
function Mfcc(mel::MelSpec; kwargs...)
    coeffs, info = _cepstrum(mel; kwargs...)
    Mfcc(coeffs, info)
end

# ---------------------------------------------------------------------------- #
#                                   get gtcc                                   #
# ---------------------------------------------------------------------------- #
"""
    Gtcc(erb::ErbSpec; ncoeffs, rect, dither) -> Gtcc

Compute Gammatone Frequency Cepstral Coefficients (GtCCs) from an ERB spectrogram.

Applies the same DCT-based cepstral pipeline as [`Mfcc`](@ref) to the output of
a gammatone filterbank ([`ErbSpec`](@ref)), yielding a biologically-motivated
alternative to standard MFCCs.

# Arguments
- `erb::ErbSpec`: Input ERB-scale gammatone spectrogram. Only [`ErbSpec`](@ref)
  is accepted; passing a [`MelSpec`](@ref) or [`BarkSpec`](@ref) will throw an
  `ArgumentError`.

# Keywords
- `ncoeffs::Int64`: Number of cepstral coefficients to retain
  (default: `nbands ÷ 2`).
- `rect::Base.Callable=mlog`: Rectification applied before DCT
  (`mlog` or `cubic_root`).
- `dither::Bool=false`: Clamp near-zero values before rectification.

# Examples
```julia
frames = Frames(audio; winsize=512, winstep=256)
stft   = Stft(frames; spectrum=power)
erb    = ErbSpec(stft; nbands=32)

gtcc = Gtcc(erb; ncoeffs=13)
```

# See also
[`Mfcc`](@ref), [`ErbSpec`](@ref), [`Delta`](@ref)
"""
function Gtcc(erb::ErbSpec; kwargs...)
    coeffs, info = _cepstrum(erb; kwargs...)
    Gtcc(coeffs, info)
end

# guard: reject non-ErbSpec input at compile-friendly runtime
function Gtcc(spec::AbstractSpectrogram; kwargs...)
    throw(ArgumentError(
        "Gtcc requires an ErbSpec (gammatone filterbank) input, got $(typeof(spec)). " *
        "Use Mfcc for MelSpec or BarkSpec inputs."
    ))
end