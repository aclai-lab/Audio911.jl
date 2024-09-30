# ---------------------------------------------------------------------------- #
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
struct MelSpecSetup
    nbands::Int
    dbscale::Bool
end

struct MelSpecData{T<:AbstractFloat}
    spec::AbstractArray{T}
    freq::AbstractVector{T}
end

struct MelSpec
    sr::Int
    setup::MelSpecSetup
    data::MelSpecData
end

function _get_melspec(
        x::AbstractArray{T},
        sr::Int;
        fbank::AbstractArray{T},
        x_length::Int,
        win_length::Int,
        noverlap::Int,
        freq::AbstractVector{T},
        nbands::Int,
        dbscale::Bool = false
) where T <: AbstractFloat
    hop_length = win_length - noverlap
    num_hops = floor(Int, (x_length - win_length) / hop_length) + 1

    # apply filterbank
    spec = reshape(fbank * x, nbands, num_hops)

    # scale to dB if requested
    if dbscale
        spec = log10.(no_zero.(spec))
    end

    MelSpec(sr, MelSpecSetup(nbands, dbscale), MelSpecData(spec, freq))
end

function Base.show(io::IO, mel_spec::MelSpec)
    print(io, "MelSpec(")
    print(io, "sr=$(mel_spec.sr), ")
    print(io, "db_scale=$(mel_spec.setup.dbscale), ")
    print(io, "spec=$(size(mel_spec.data.spec)), ")
    print(io, "freq=$(length(mel_spec.data.freq)) frequencies)")
end

function Base.display(mel_spec::MelSpec)
    heatmap(1:size(mel_spec.data.spec, 2), 1:size(mel_spec.data.spec, 1), mel_spec.data.spec;
        ylabel="Bands",
        xlabel="Frame",
        title="Mel Spectrogram",
    )
end

function get_melspec(source::Stft; fbank::MelFb, kwargs...)
    _get_melspec(
        source.data.spec,
        source.sr;
        fbank=fbank.data.fbank,
        x_length=source.x_length,
        win_length=source.setup.nwin, 
        noverlap=source.setup.noverlap,
        freq=fbank.data.freq,
        nbands=fbank.setup.nbands,
        kwargs...)
end