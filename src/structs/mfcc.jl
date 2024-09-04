# ---------------------------------------------------------------------------- #
#                                     mfcc                                     #
# ---------------------------------------------------------------------------- #
struct MfccSetup
    ncoeffs::Int64
	rectification::Symbol
	dither::Bool
end

struct MfccData
    spec::AbstractArray{<:AbstractFloat}
end

struct Mfcc
    sr::Int64
    setup::MfccSetup
    data::MfccData
end

function create_DCT_matrix(mel_coeffs::Int64)
    # create DCT matrix
    matrix = zeros(Float64, mel_coeffs, mel_coeffs)
    s0 = sqrt(1 / mel_coeffs)
    s1 = sqrt(2 / mel_coeffs)
    piCCast = 2 * pi / (2 * mel_coeffs)

    matrix[1, :] .= s0
    for k in 1:mel_coeffs, n in 2:mel_coeffs
        matrix[n, k] = s1 * cos(piCCast * (n - 1) * (k - 0.5))
    end

    matrix
end

function _get_mfcc(;
    spec::AbstractArray{Float64},
    sr::Int64,
    nbands::Int64,
    ncoeffs::Int64 = round(Int, nbands / 2),
    rect::Symbol = :log, # :log, :cubic_root
    dither::Bool = false,
)
    dither ? spec[spec .< 1e-8] .= 1e-8 : spec[spec .== 0] .= floatmin(Float64)

    # Design DCT matrix
    DCTmatrix = create_DCT_matrix(nbands)

    # apply DCT matrix
    if (rect == :log)
        coeffs = DCTmatrix * log10.(spec)
    elseif (rect == :cubic_root)
        # apply DCT matrix
        coeffs = DCTmatrix * spec .^ (1 / 3)
    else
        @warn("Unknown $rect DCT matrix rectification, defaulting to log.")
        coeffs = DCTmatrix * log10.(spec)
    end

    Mfcc(sr, MfccSetup(ncoeffs, rect, dither), MfccData(coeffs[1:ncoeffs, :]))
end

function Base.show(io::IO, mfcc::Mfcc)
    println(io, "Mfcc:")
    println(io, "  Sample rate: $(mfcc.sr) Hz")
    println(io, "  Number of coefficients: $(mfcc.setup.nbands)")
    println(io, "  Dither: $(mfcc.setup.dither)")
    println(io, "  Rectification: $(mfcc.setup.rectification)")
    println(io, "  MFCC shape: $(size(mfcc.data.spec))")
end

function Base.display(mfcc::Mfcc)    
    heatmap(1:size(mfcc.data.spec, 2), 1:size(mfcc.data.spec, 1), mfcc.data.spec;
        ylabel="Bands",
        xlabel="Frame",
        title="Mfcc",
    )
end

function get_mfcc(;source = nothing, kwargs...)
    if source isa MelSpec
        _get_mfcc(; spec=source.data.spec, sr=source.sr, nbands=source.setup.nbands, kwargs...)
    elseif source isa Cwt
        _get_mfcc(; spec=source.data.spec, sr=source.sr, nbands=length(source.setup.freq), kwargs...)
    else
        error("Unsupported type for spec: typeof(source) = $(typeof(source))")
    end
end