# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
struct StftSetup
    nfft::Int
    wintype::Tuple{Symbol, Symbol}
    nwin::Int
    noverlap::Int
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

const NORM_FUNCS = Dict(:power => x -> (@. real((x * conj(x)))), :magnitude => x -> (@. abs(x)))

function _get_stft(
    x::AbstractVector{T},
    sr::Int;
    nfft::Union{Int, AbstractFloat, Nothing} = nothing,
    nwin::Union{Int, AbstractFloat, Nothing} = nothing,
    noverlap::Union{Int, AbstractFloat, Nothing} = nothing,
    wintype::Tuple{Symbol, Symbol} = (:hann, :periodic),
    norm::Symbol = :power, # :none, :power, :magnitude
) where T<:AbstractFloat
    typeof(nfft) <: AbstractFloat && begin nfft = round(Int, nfft * sr) end
    typeof(nwin) <: AbstractFloat && begin nwin = round(Int, nwin * sr) end
    typeof(noverlap) <: AbstractFloat && begin noverlap = round(Int, noverlap * sr) end

    # apply default parameters if not provided
    nfft = nfft !== nothing ? nfft : sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512
    nwin = nwin !== nothing ? nwin : nfft
    noverlap = noverlap !== nothing ? noverlap : round(Int, nwin * 0.5)

    @assert noverlap < nwin "Overlap length must be smaller than nwin: $nwin."
    @assert nwin ≤ nfft "FFT window size smaller than actual window size is highly discuraged."
    @assert haskey(NORM_FUNCS, norm) "Unknown spectrum_type: $norm."

    frames, win = _get_frames(x, wintype, nwin, noverlap)    
    wframes = nwin < nfft ? vcat(frames, zeros(Float64, nfft - nwin, size(frames, 2))) .* win : frames .* win

    nout = (nfft >> 1)+1
    spec = zeros(T, nout, size(wframes, 2))
    tmp = Vector{ComplexF64}(undef, nout)

    plan = plan_rfft(1:nfft)
    offset = 0

    @inbounds @simd for i in eachcol(wframes)
        mul!(tmp, plan, i)
        copyto!(spec, offset+1, NORM_FUNCS[norm](tmp), 1, nout)
        offset += nout
    end
    
    freq = sr / nfft * (StepRangeLen(1:nout) .- 1)

    Stft(sr, size(x, 1), StftSetup(nfft, wintype, nwin, noverlap, norm), StftData(spec, freq, win, frames))
end

function Base.show(io::IO, stft::Stft)
    print(io, "Stft(")
    print(io, "nfft=$(stft.setup.nfft), ")
    print(io, "win_type=$(stft.setup.wintype), ")
    print(io, "win_length=$(stft.setup.nwin), ")
    print(io, "overlap_length=$(stft.setup.noverlap), ")
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

get_stft(source::Audio; kwargs...) = _get_stft(source.data, source.sr; kwargs...)