# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct MfccSetup <: AbstractSetup
    sr      :: Int64
    ncoeffs :: Int64
	rect    :: Base.Callable
	dither  :: Bool
end

# ---------------------------------------------------------------------------- #
#                                     mfcc                                     #
# ---------------------------------------------------------------------------- #
struct Mfcc{T} <: AbstractSpectrogram
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
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    Base.eltype(::Mfcc{T}) -> Type

Return the element type of the mfcc spectrogram data.
"""
Base.eltype(::Mfcc{T}) where {T} = T

"""
    get_data(m::Mfcc) -> Matrix

Get the mfcc spectrogram data matrix, transposed to (nframes × nbands).
"""
@inline get_data(m::Mfcc)  = m.spec'

"""
    get_setup(m::Mfcc) -> MelSpecSetup

Get the configuration metadata for the mfcc spectrogram.
"""
@inline get_setup(m::Mfcc) = m.info

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
#                                   get mfcc                                   #
# ---------------------------------------------------------------------------- #
"""
    Mfcc(spec::MelSpec, sr::Int64; kwargs...) -> Mfcc

Compute Mel-Frequency Cepstral Coefficients (MFCCs) from a mel spectrogram.

MFCCs are derived by applying a Discrete Cosine Transform (DCT) to the logarithmic 
mel-frequency energies. They provide a compact representation of the spectral envelope 
and are widely used in speech recognition, speaker identification, and music analysis.

# Arguments
- `spec::MelSpec`: Input mel spectrogram
- `sr::Int64`: Sample rate in Hz

# Keyword Arguments
- `nbands::Int64`: Number of mel-frequency bands in the input spectrogram
- `ncoeffs::Int64`: Number of cepstral coefficients to retain (default: `nbands/2`)
- `rect::Base.Callable`: Rectification function applied before DCT (default: `mlog`)
  - `mlog`: Logarithmic scaling (recommended)
  - `cubic_root`: Cube root scaling (alternative)
- `dither::Bool`: Apply dithering to prevent log(0) by setting small values to 1e-8 
  (default: `false`). When `false`, zeros are replaced with `floatmin(Float64)`.

# Examples
```julia
# Standard 13 MFCC coefficients
melspec = MelSpec(stft; nbands=40)
mfcc = Mfcc(melspec, 16000; nbands=40, ncoeffs=13)

# Using cube root rectification with dithering
mfcc = Mfcc(melspec, 16000; nbands=40, ncoeffs=20, rect=cubic_root, dither=true)
```

# See also
[`MelSpec`](@ref), [`mlog`](@ref), [`cubic_root`](@ref), [`create_dct_matrix`](@ref)
"""
function Mfcc(
    spec    :: MelSpec,
    sr      :: Int64;
    nbands  :: Int64,
    ncoeffs :: Int64=round(Int, nbands / 2),
    rect    :: Base.Callable=mlog, # mlog, cubic_root
    dither  :: Bool=false,
)::Mfcc
    dither ? spec[spec .< 1e-8] .= 1e-8 : spec[spec .== 0] .= floatmin(Float64)

    dct_matrix = _create_dct_matrix(nbands)
    coeffs = rect(dct_matrix, spec)

    info = MfccSetup(sr, ncoeffs, rect, dither)

    Mfcc(coeffs[1:ncoeffs, :], info)
end
