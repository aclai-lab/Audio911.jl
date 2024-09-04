# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
struct StftSetup
    nfft::Int
    wintype::Tuple{Symbol, Symbol}
    winlength::Int
    overlaplength::Int
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
        nfft::Union{Int, Time, Nothing} = nothing,
        winlength::Union{Int, Time, Nothing} = nothing,
        overlaplength::Union{Int, Time, Nothing} = nothing,
        wintype::Tuple{Symbol, Symbol} = (:hann, :periodic),
        norm::Symbol = :power, # :none, :power, :magnitude, :pow2mag
)
    # ms to sample conversion
    typeof(nfft) <: Time && begin nfft = round(Int, ustrip(Int64, u"ms", nfft) * sr * 0.001) end
    typeof(winlength) <: Time && begin winlength = round(Int, ustrip(Int64, u"ms", winlength) * sr * 0.001) end
    typeof(overlaplength) <: Time && begin overlaplength = round(Int, ustrip(Int64, u"ms", overlaplength) * sr * 0.001) end

    # apply default parameters if not provided
    nfft = nfft !== nothing ? nfft : sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512
    winlength = winlength !== nothing ? winlength : nfft
    overlaplength = overlaplength !== nothing ? overlaplength : round(Int, nfft / 2)

    @assert overlaplength<nfft "Overlap length must be smaller than nfft."

    x_length = size(x, 1)
    frames, win, wframes, _, _ = _get_frames(x, wintype, winlength, overlaplength)

    if winlength < nfft
        wframes = vcat(wframes, zeros(Float64, nfft - winlength, size(wframes, 2)))
    elseif winlength > nfft
        @error("FFT window size smaller than actual window size is highly discuraged.")
    end

    # take one side
    if mod(nfft, 2) == 0
        oneside = 1:Int(nfft / 2 + 1)   # even
    else
        oneside = 1:Int((nfft + 1) / 2)  # odd
    end

    # normalize
    norm_funcs = Dict(
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => x -> sqrt.(real.((x .* conj.(x))))
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, norm) "Unknown spectrum_type: $norm."

    spec = norm_funcs[norm](fft(wframes, (1,))[oneside, :])
    freq = (sr / nfft) * (oneside .- 1)

    Stft(sr, x_length, StftSetup(nfft, wintype, winlength, overlaplength, norm), StftData(spec, freq, win, frames))
end

function Base.show(io::IO, stft::Stft)
    print(io, "Stft(")
    print(io, "nfft=$(stft.setup.nfft), ")
    print(io, "win_type=$(stft.setup.wintype), ")
    print(io, "win_length=$(stft.setup.winlength), ")
    print(io, "overlap_length=$(stft.setup.overlaplength), ")
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

get_stft(; source::Audio, kwargs...) = _get_stft(; x = source.data, sr = source.sr, kwargs...)