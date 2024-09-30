# ---------------------------------------------------------------------------- #
#                                    deltas                                    #
# ---------------------------------------------------------------------------- #
struct DeltasSetup
    dlength::Int
    transpose::Bool
end

struct DeltasData{T<:AbstractFloat}
    dspec::AbstractArray{T}
    ddspec::AbstractArray{T}
end

struct Deltas
    sr::Int
    setup::DeltasSetup
    data::DeltasData
end

function audiodelta(x::AbstractArray{<:AbstractFloat}, dlength::Int, transpose::Bool)
    # define window shape
    m = floor(Int, dlength / 2)
    b = collect(m:-1:(-m)) ./ sum((1:m) .^ 2)

    transpose ? filt(b, 1.0, x')' : filt(b, 1.0, x)
end

function _get_deltas(
    source::AbstractArray{<:AbstractFloat}, 
    sr::Int;
    dlength::Int = 9, 
    transpose::Bool = false
)
    delta = audiodelta(source, dlength, transpose)
    ddelta = audiodelta(delta, dlength, transpose)

    Deltas(sr, DeltasSetup(dlength, transpose), DeltasData(delta, ddelta))
end

# function Base.show(io::IO, deltas::Deltas)
#     println(io, "Deltas:")
#     println(io, "  Delta window length: $(deltas.d_length)")
#     println(io, "  Delta matrix type: $(deltas.d_matrix)")
#     println(io, "  Delta shape: $(size(deltas.delta))")
#     println(io, "  Delta-delta shape: $(size(deltas.ddelta))")
#     println(io, "  Frequency range: $(first(deltas.freq)) - $(last(deltas.freq)) Hz")
# end

function Base.display(deltas::Deltas)    
    # Plot delta features
    d1 = heatmap(1:size(deltas.data.dspec, 2), 1:size(deltas.data.dspec, 1), deltas.data.dspec;
        ylabel="Bands",
        xlabel="Frame",
        title="Delta",
    )
    
    # Plot delta-delta features
    d2 = heatmap(1:size(deltas.data.ddspec, 2), 1:size(deltas.data.ddspec, 1), deltas.data.ddspec;
        ylabel="Bands",
        xlabel="Frame",
        title="Delta-Delta",
    )

    plot(d1, d2, layout=(2,1), size=(800,1000), link=:x)
    display(current())
end

function get_deltas(source::Mfcc; kwargs...)
     _get_deltas(source.data.spec, source.sr; kwargs...)
end
