# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
struct StftSetup
    nfft::Int
    win_type::Tuple{Symbol, Symbol}
    win_length::Int
    overlap_length::Int
    norm::Symbol 
end

struct StftData
    spec::AbstractArray{<:AbstractFloat}
    freq::StepRangeLen{<:AbstractFloat}
    win::AbstractVector{<:AbstractFloat}
    frames::AbstractArray{<:AbstractFloat}
end

struct Stft
    sr::Int64
    x_length::Int64
    setup::StftSetup
    data::StftData
end

function _get_stft(;
        x::AbstractVector{<:AbstractFloat},
        sr::Int,
        # nfft::Int = sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512,
        # win_length::Int = nfft,
        # overlap_length::Int = round(Int, nfft / 2),
        nfft::Union{Int, Time, Nothing} = nothing,
        win_length::Union{Int, Time, Nothing} = nothing,
        overlap_length::Union{Int, Time, Nothing} = nothing,
        win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
        norm::Symbol = :power, # :none, :power, :magnitude, :pow2mag
)
    # ms to sample conversion
    typeof(nfft) <: Time && begin nfft = round(Int, ustrip(Int64, u"ms", nfft) * sr * 0.001) end
    typeof(win_length) <: Time && begin win_length = round(Int, ustrip(Int64, u"ms", win_length) * sr * 0.001) end
    typeof(overlap_length) <: Time && begin overlap_length = round(Int, ustrip(Int64, u"ms", overlap_length) * sr * 0.001) end

    # apply default parameters if not provided
    nfft = nfft !== nothing ? nfft : sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512
    win_length = win_length !== nothing ? win_length : nfft
    overlap_length = overlap_length !== nothing ? overlap_length : round(Int, nfft / 2)

    @assert overlap_length<nfft "Overlap length must be smaller than nfft."

    x_length = size(x, 1)
    frames, win, wframes, _, _ = _get_frames(x, win_type, win_length, overlap_length)

    if win_length < nfft
        wframes = vcat(wframes, zeros(Float64, nfft - win_length, size(wframes, 2)))
    elseif win_length > nfft
        @error("FFT window size smaller than actual window size is highly discuraged.")
    end

    # take one side
    if mod(nfft, 2) == 0
        one_side = 1:Int(nfft / 2 + 1)   # even
    else
        one_side = 1:Int((nfft + 1) / 2)  # odd
    end

    # normalize
    norm_funcs = Dict(
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => x -> sqrt.(real.((x .* conj.(x))))
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, norm) "Unknown spectrum_type: $norm."

    spec = norm_funcs[norm](fft(wframes, (1,))[one_side, :])
    freq = (sr / nfft) * (one_side .- 1)

    Stft(sr, x_length, StftSetup(nfft, win_type, win_length, overlap_length, norm), StftData(spec, freq, win, frames))
end

function Base.show(io::IO, stft::Stft)
    print(io, "Stft(")
    print(io, "nfft=$(stft.setup.nfft), ")
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
    _get_stft(; x = audio.data, sr = audio.sr, kwargs...)
end