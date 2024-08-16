# ---------------------------------------------------------------------------- #
#                            linear spectrogram                                #
# ---------------------------------------------------------------------------- #
mutable struct LinSpec
	# setup
    sr::Int64
    freq_range::Tuple{Int64, Int64}
	db_scale::Bool
	# data
	spec::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function LinSpec(;
    sr,
    freq_range = (0, round(Int, sr / 2)),
	db_scale = false,
	spec = nothing,
	freq = nothing,
)
    lspec = LinSpec(
        sr,
        freq_range,
        db_scale,
        spec,
        freq
    )
	return 	lspec
end

# helper function to avoid log of zero
no_zero = x -> x == 0 ? floatmin(Float64) : x

function _get_lin_spec(;
    x::AbstractArray{Float64},
    xfreq::AbstractVector{Float64},
    lspec::LinSpec
)
    # trim to desired frequency range
    x_range = findall(lspec.freq_range[1] .<= xfreq .<= lspec.freq_range[2])
    lspec.spec, lspec.freq = x[x_range, :], xfreq[x_range]

    # scale to dB if requested
    if lspec.db_scale
        lspec.spec = log10.(no_zero.(lspec.spec))
    end

    return lspec
end

function Base.show(io::IO, lspec::LinSpec)
    print(io, "LinSpec(")
    print(io, "freq_range=$(lspec.freq_range), ")
    print(io, "db_scale=$(lspec.db_scale), ")
    if isnothing(lspec.spec)
        print(io, "spec=nothing, ")
    else
        print(io, "spec=$(size(lspec.spec)), ")
    end
    if isnothing(lspec.freq)
        print(io, "freq=nothing)")
    else
        print(io, "freq=$(length(lspec.freq)) frequencies)")
    end
end

function Base.display(lspec::LinSpec)
    isnothing(lspec.spec) && error("Linear spectrogram has not been computed yet.")
    
    # define frequency ticks
    min_freq, max_freq = extrema(lspec.freq)
    freq_ticks = 10 .^ range(log10(min_freq + eps()), log10(max_freq + eps()), length=5)
    freq_ticks = round.(freq_ticks, digits=1)
    
    heatmap(1:size(lspec.spec, 2), lspec.freq, lspec.spec;
        ylabel="Frequency (Hz)",
        xlabel="Frame",
        title="CWT Spectrogram",
        yaxis=:log,
        yticks=(freq_ticks, string.(freq_ticks))
    )
end

function get_linspec(;
    source = nothing,
    kwargs...
)
    if source isa Stft
        _get_lin_spec(x=source.data.spec, xfreq=source.data.freq, lspec=LinSpec(;source.sr, kwargs...))
    elseif source isa Cwt
        _get_lin_spec(x=source.data.spec, xfreq=source.data.freq, lspec=LinSpec(;source.sr, freq_range=round.(Int, extrema(source.freq)), kwargs...))
    else
        error("Unsupported type for spec: typeof(source) = $(typeof(source))")
    end
end
