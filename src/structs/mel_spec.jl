# ---------------------------------------------------------------------------- #
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
mutable struct MelSpec
	# setup
    sr::Int64
    nbands::Int64
    db_scale::Bool
	# data
	spec::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function MelSpec(;
    sr,
    nbands,
    db_scale = false,
    spec = nothing,
    freq = nothing,
)
    mel_spec = MelSpec(
        sr,
        nbands,
        db_scale,
        spec,
        freq
    )
	return 	mel_spec
end

function _get_melspec(;
        x::AbstractArray{Float64},
        fbank::AbstractArray{Float64},
        x_length::Int64,
        win_length::Int64,
        overlap_length::Int64,
        freq::AbstractVector{Float64},
        mel_spec::MelSpec
)
    hop_length = win_length - overlap_length
    num_hops = floor(Int, (x_length - win_length) / hop_length) + 1

    # apply filterbank
    mel_spec.spec = reshape(fbank * x, mel_spec.nbands, num_hops)
    mel_spec.freq = freq

    # scale to dB if requested
    if mel_spec.db_scale
        mel_spec.spec = log10.(no_zero.(mel_spec.spec))
    end

    return mel_spec
end # melSpectrogram

function Base.show(io::IO, mel_spec::MelSpec)
    print(io, "MelSpec(")
    print(io, "sr=$(mel_spec.sr), ")
    print(io, "db_scale=$(mel_spec.db_scale), ")
    if isnothing(mel_spec.spec)
        print(io, "spec=nothing, ")
    else
        print(io, "spec=$(size(mel_spec.spec)), ")
    end
    if isnothing(mel_spec.freq)
        print(io, "freq=nothing)")
    else
        print(io, "freq=$(length(mel_spec.freq)) frequencies)")
    end
end

function Base.display(mel_spec::MelSpec)
    isnothing(mel_spec.spec) && error("Mel spectrogram has not been computed yet.")
    
    heatmap(1:size(mel_spec.spec, 2), 1:size(mel_spec.spec, 1), mel_spec.spec;
        ylabel="Bands",
        xlabel="Frame",
        title="Mel Spectrogram",
    )
end

function get_melspec(;
    stft::Stft,
    fbank::MelFbank,
    kwargs...
)
    _get_melspec(;
        x=stft.spec,
        fbank=fbank.fbank,
        x_length=stft.x_length,
        win_length=stft.win_length, 
        overlap_length=stft.overlap_length,
        freq=fbank.freq,
        mel_spec=MelSpec(; sr=stft.sr, nbands=fbank.nbands, kwargs...)
    )
end