# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
struct StftSetup
    stft_length::Int
    win_type::Tuple{Symbol, Symbol}
    win_length::Int
    overlap_length::Int
    norm::Symbol 
end

struct StftData
    spec::AbstractArray{<:AbstractFloat}
    freq::AbstractVector{<:AbstractFloat}
    win::AbstractVector{<:AbstractFloat}
    frames::AbstractArray{<:AbstractFloat}
end

struct Stft
    # setup
    sr::Int64
    x_length::Int64
    setup::StftSetup
    data::StftData
end

function _get_stft(;
        x::AbstractVector{<:AbstractFloat},
        sr::Int,
        x_length::Int,
        stft_length::Int = sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512,
        win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
        win_length::Int = stft_length,
        overlap_length::Int = round(Int, stft_length / 2),
        stft_norm::Symbol = :power, # :none, :power, :magnitude, :pow2mag
        _whatever...
)
    frames, win, wframes, _, _ = _get_frames(x, win_type, win_length, overlap_length)

    if win_length < stft_length
        wframes = vcat(wframes, zeros(Float64, stft_length - win_length, size(wframes, 2)))
    elseif win_length > stft_length
        @error("FFT window size smaller than actual window size is highly discuraged.")
    end

    # take one side
    if mod(stft_length, 2) == 0
        one_side = 1:Int(stft_length / 2 + 1)   # even
    else
        one_side = 1:Int((stft_length + 1) / 2)  # odd
    end

    # normalize
    norm_funcs = Dict(
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => x -> sqrt.(real.((x .* conj.(x))))
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, stft_norm) "Unknown spectrum_type: $stft_norm."

    spec = norm_funcs[stft_norm](fft(wframes, (1,))[one_side, :])
    freq = (sr / stft_length) * (one_side .- 1)

    Stft(sr, x_length, StftSetup(stft_length, win_type, win_length, overlap_length, stft_norm), StftData(spec, freq, win, frames))
end

function Base.show(io::IO, stft::Stft)
    print(io, "Stft(")
    print(io, "stft_length=$(stft.setup.stft_length), ")
    print(io, "win_type=$(stft.setup.win_type), ")
    print(io, "win_length=$(stft.setup.win_length), ")
    print(io, "overlap_length=$(stft.setup.overlap_length), ")
    print(io, "norm=:$(stft.setup.norm), ")
    print(io, "spec=$(size(stft.data.spec)))")
end

function Plots.plot(stft::Stft)
    heatmap(1:size(stft.data.spec, 2), stft.data.freq, stft.data.spec;
        ylabel = "Frequency (Hz)",
        xlabel = "Frame",
        title = "STFT Spectrogram",
        color = :viridis,
        clim = (0, 1)
    )
end

function get_stft(;
        audio::Audio,
        kwargs...
)
    _get_stft(; x = audio.data, sr = audio.sr, x_length = size(audio.data, 1), kwargs...)
end