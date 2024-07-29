# ---------------------------------------------------------------------------- #
#                                     Mfcc                                     #
# ---------------------------------------------------------------------------- #
mutable struct Mfcc
	# setup
    sr::Int64
    ncoeffs::Int64
	rectification::Symbol
	dither::Bool
	# data
	mfcc::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function Mfcc(;
	sr,
    nbands,
    ncoeffs = nothing,
    rectification = :log, # :log, :cubic_root
    dither = false,
	mfcc = nothing,
	freq = nothing,
)
    if isnothing(ncoeffs)
        ncoeffs = round(Int, nbands / 2)
    end

	mfcc = Mfcc(
        sr,
        ncoeffs,
		rectification,
		dither,
		mfcc,
		freq
	)
	return mfcc
end

mutable struct Deltas
	# setup
    sr::Int64
    d_length::Int64
    d_matrix::Symbol
	# data
    delta::Union{Nothing, AbstractArray{Float64}}
    ddelta::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function Deltas(;
	sr,
    d_length = 9,
    d_matrix = :standard, # :standard, :transposed
	delta = nothing,
	ddelta = nothing,
	freq = nothing,
)
	deltas = Deltas(
        sr,
        d_length,
		d_matrix,
		delta,
		ddelta,
		freq
	)
	return deltas
end

# ---------------------------------------------------------------------------- #
#                                     Dct II                                   #
# ---------------------------------------------------------------------------- #
function create_DCT_matrix(
        mel_coeffs::Int64,
)
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

# ---------------------------------------------------------------------------- #
#                                     Deltas                                   #
# ---------------------------------------------------------------------------- #
function audioDelta(
        x::AbstractArray{Float64},
        win_length::Int64,
        source::Symbol = :standard
)

    # define window shape
    m = Int(floor(win_length / 2))
    b = collect(m:-1:(-m)) ./ sum((1:m) .^ 2)

    if source == :transposed
        filt(b, 1.0, x')'   #:audioflux setting
    else
        filt(b, 1.0, x)     #:matlab setting
    end
end

# ---------------------------------------------------------------------------- #
#                                     Mfcc                                     #
# ---------------------------------------------------------------------------- #
function _get_mfcc(;
        spec::AbstractArray{Float64},
        nbands::Int64,
        freq::AbstractVector{Float64},
        mfcc::Mfcc,
)
    mfcc.dither ? spec[spec .< 1e-8] .= 1e-8 : spec[spec .== 0] .= floatmin(Float64)

    # Design DCT matrix
    DCTmatrix = create_DCT_matrix(nbands)

    # apply DCT matrix
    if (mfcc.rectification == :log)
        coeffs = DCTmatrix * log10.(spec)
    elseif (mfcc.rectification == :cubic_root)
        # apply DCT matrix
        coeffs = DCTmatrix * spec .^ (1 / 3)
    else
        @warn("Unknown $mfcc.rectification DCT matrix rectification, defaulting to log.")
        coeffs = DCTmatrix * log10.(spec)
    end

    # reduce to mfcc coefficients
    mfcc.mfcc = coeffs[1:mfcc.ncoeffs, :]
    mfcc.freq = freq

    return mfcc
end

function Base.show(io::IO, mfcc::Mfcc)
    println(io, "Mfcc:")
    println(io, "  Sample rate: $(mfcc.sr) Hz")
    println(io, "  Number of coefficients: $(mfcc.ncoeffs)")
    println(io, "  Dither: $(mfcc.dither)")
    println(io, "  Rectification: $(mfcc.rectification)")
    println(io, "  MFCC shape: $(size(mfcc.mfcc))")
    println(io, "  Frequency range: $(first(mfcc.freq)) - $(last(mfcc.freq)) Hz")
end

function Base.display(mfcc::Mfcc)
    isnothing(mfcc.mfcc) && error("Mfcc has not been computed yet.")
    
    heatmap(1:size(mfcc.mfcc, 2), 1:size(mfcc.mfcc, 1), mfcc.mfcc;
        ylabel="Bands",
        xlabel="Frame",
        title="Mfcc",
    )
end

function get_mfcc(;
    source = nothing,
    kwargs...
)
    if source isa MelSpec
        _get_mfcc(spec=source.spec, nbands=source.nbands, freq=source.freq, mfcc=Mfcc(; sr=source.sr, nbands=source.nbands, kwargs...))
    elseif source isa Cwt
        _get_mfcc(spec=source.spec, nbands=length(source.freq), freq=source.freq, mfcc=Mfcc(; sr=source.sr, nbands=length(source.freq), kwargs...))
    else
        error("Unsupported type for spec: typeof(source) = $(typeof(source))")
    end
end

function _get_deltas(;
        source::AbstractArray{Float64},
        freq::AbstractVector{Float64},
        deltas::Deltas,
)
    deltas.delta = audioDelta(
        source, deltas.d_length, deltas.d_matrix)
    deltas.ddelta = audioDelta(
        deltas.delta, deltas.d_length, deltas.d_matrix)

    deltas.freq = freq

    return deltas
end

function Base.show(io::IO, deltas::Deltas)
    println(io, "Deltas:")
    println(io, "  Delta window length: $(deltas.d_length)")
    println(io, "  Delta matrix type: $(deltas.d_matrix)")
    println(io, "  Delta shape: $(size(deltas.delta))")
    println(io, "  Delta-delta shape: $(size(deltas.ddelta))")
    println(io, "  Frequency range: $(first(deltas.freq)) - $(last(deltas.freq)) Hz")
end

function Base.display(deltas::Deltas)
    isnothing(deltas.delta) && error("Delta has not been computed yet.")
    isnothing(deltas.ddelta) && error("Delta-Delta has not been computed yet.")
    
    # Plot delta features
    d1 = heatmap(1:size(deltas.delta, 2), 1:size(deltas.delta, 1), deltas.delta;
        ylabel="Bands",
        xlabel="Frame",
        title="Delta",
    )
    
    # Plot delta-delta features
    d2 = heatmap(1:size(deltas.ddelta, 2), 1:size(deltas.ddelta, 1), deltas.ddelta;
        ylabel="Bands",
        xlabel="Frame",
        title="Delta-Delta",
    )

    plot(d1, d2, layout=(2,1), size=(800,1000), link=:x)
    display(current())
end

function get_deltas(;
    source::Mfcc,
    kwargs...
)
    _get_deltas(source=source.mfcc, freq=source.freq, deltas=Deltas(; sr=source.sr, kwargs...))
end
