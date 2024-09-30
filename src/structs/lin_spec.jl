# ---------------------------------------------------------------------------- #
#                            linear spectrogram                                #
# ---------------------------------------------------------------------------- #
struct LinSpecSetup
    freqrange::Tuple{Int, Int}
    dbscale::Bool
end

struct LinSpecData{T<:AbstractFloat}
    spec::AbstractArray{<:T}
    freq::AbstractVector{<:T}
end

struct LinSpec
    sr::Int
    setup::LinSpecSetup
    data::LinSpecData
end

# helper function to avoid log of zero
no_zero = x -> x == 0 ? floatmin(Float64) : x

function _get_lin_spec(
    x::AbstractArray{T},
    sr::Int;
    xfreq::AbstractVector{T},
    freqrange::Tuple{Int, Int} = (0, round(Int, sr / 2)),
	dbscale::Bool = false,
) where {T<:AbstractFloat}
    # trim to desired frequency range
    x_range = findall(freqrange[1] .<= xfreq .<= freqrange[2])
    spec, freq = x[x_range, :], xfreq[x_range]

    # scale to dB if requested
    if dbscale
        spec = log10.(no_zero.(spec))
    end

    LinSpec(sr, LinSpecSetup(freqrange, dbscale), LinSpecData(spec, freq))
end

function Base.show(io::IO, lspec::LinSpec)
    print(io, "LinSpec(")
    print(io, "freq_range=$(lspec.setup.freqrange), ")
    print(io, "db_scale=$(lspec.setup.dbscale), ")
end

function Base.display(lspec::LinSpec)    
    # define frequency ticks
    min_freq, max_freq = extrema(lspec.data.freq)
    freq_ticks = 10 .^ range(log10(min_freq + eps()), log10(max_freq + eps()), length=5)
    freq_ticks = round.(freq_ticks, digits=1)
    
    heatmap(1:size(lspec.data.spec, 2), lspec.data.freq, lspec.data.spec;
        ylabel="Frequency (Hz)",
        xlabel="Frame",
        title="Linear Spectrogram",
        yaxis=:log,
        yticks=(freq_ticks, string.(freq_ticks))
    )
end

function get_linspec(source = nothing; kwargs...)
    if source isa Stft
        _get_lin_spec(source.data.spec, source.sr; xfreq=source.data.freq, kwargs...)
    elseif source isa Cwt
        _get_lin_spec(source.data.spec, source.sr; xfreq=source.data.freq, kwargs...)
        ### BROKEN >>> lspec=LinSpec(;source.sr, freq_range=round.(Int, extrema(source.freq)), kwargs...))
    else
        error("Unsupported type for spec: typeof(source) = $(typeof(source))")
    end
end
