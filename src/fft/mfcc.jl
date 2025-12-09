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
    mlog -> Function

Create a function that applies base-10 logarithmic scaling to spectral data.

Returns a function that performs element-wise logarithmic transformation on a matrix,
commonly used as a rectification method in MFCC computation to compress the dynamic
range of mel-frequency energies.
"""
mlog(x, y) =  x * log10.(y)

"""
    cubic_root -> Function

Create a function that applies cube root scaling to spectral data.

Returns a function that performs element-wise cube root transformation on a matrix,
used as an alternative rectification method in MFCC computation. Provides compression
similar to logarithmic scaling but with different mathematical properties.
"""
cubic_root(x, y) = x * y .^ (1 / 3)


function create_dct_matrix(mel_coeffs::Int64)
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
function Mfcc(
    spec    :: MelSpec,
    sr      :: Int64;
    nbands  :: Int64,
    ncoeffs :: Int64=round(Int, nbands / 2),
    rect    :: Base.Callable=mlog, # mlog, cubic_root
    dither  :: Bool=false,
)::Mfcc
    dither ? spec[spec .< 1e-8] .= 1e-8 : spec[spec .== 0] .= floatmin(Float64)

    dct_matrix = create_dct_matrix(nbands)
    coeffs = rect(dct_matrix, spec)

    info = MfccSetup(sr, ncoeffs, rect, dither)

    Mfcc(coeffs[1:ncoeffs, :], info)
end


# function get_mfcc(source::Union{MelSpec, Cwt, Nothing}=nothing; kwargs...)
#     if source isa MelSpec
#         _get_mfcc(source.data.spec, source.sr; nbands=source.setup.nbands, kwargs...)
#     elseif source isa Cwt
#         _get_mfcc(source.data.spec, source.sr; nbands=length(source.setup.freq), kwargs...)
#     end
# end