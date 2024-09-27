# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
struct StftSetup
    nfft::Int
    nwin::Int
    noverlap::Int
    wintype::Tuple{Symbol, Symbol}
    norm::Symbol 
end

struct StftData{T<:AbstractFloat}
    spec::AbstractArray{T}
    freq::StepRangeLen{T}
    win::AbstractVector{T}
    frames::Base.Generator
end

struct Stft
    sr::Int64
    x_length::Int64
    setup::StftSetup
    data::StftData
end

global const NORM_FUNCS = Dict{Symbol, Function}(
    :power => (x, y) -> (@. real(x * conj(x))),
    :magnitude => (x, y) -> (@. abs(x)),
    :winpower => (x, y) -> real.(x .* conj(x)) ./ (0.5 * sum(y)^2),
    :winmagnitude => (x, y) -> abs.(x) ./ (0.5 * sum(y))
)

@kwdef struct WindowFunctions
    hann::Function = x -> 0.5 * (1 + cospi(2x))
    hamming::Function = x -> 0.54 - 0.46 * cospi(2x)
    rect::Function = x -> 1.0
end

function _get_window(wintype::Symbol, nwin::Int, winperiod::Bool)
    nwin == 1 && return [1.0]
    winfunc = getproperty(WindowFunctions(), wintype)
    winperiod && return collect(winfunc(x) for x in range(-0.5, stop=0.5, length=nwin+1)[1:end-1])
    collect(winfunc(x) for x in range(-0.5, stop=0.5, length=nwin))
end

function _get_frames(x::AbstractVector{<:AbstractFloat}, nwin::Int, noverlap::Int)
    nhop = nwin - noverlap
    nhops = div(length(x) - nwin, nhop) + 1
    @views (view(x, i:i+nwin-1) for i in 1:nhop:(nhops-1)*nhop+1)
end

_get_wframes(frames::Base.Generator, win::AbstractVector{<:AbstractFloat}) = @inbounds collect(i .* win for i in frames)

@inline function _get_stft(
    x::AbstractVector{T},
    sr::Int;
    nfft::Int = prevpow(2, sr ÷ 30),
    nwin::Int = nfft,
    noverlap::Int = round(Int, nwin * 0.5),
    wintype::Tuple{Symbol, Symbol} = (:hann, :periodic),
    norm::Symbol = :power, # :none, :power, :magnitude
    halve::Bool = true,
) where T<:AbstractFloat
    (0 ≤ noverlap ≤ nwin) || throw(DomainError((; noverlap, nwin), "Overlap length must be smaller than nwin: $nwin."))
    nfft >= nwin || throw(DomainError((; nfft, nwin), "nfft must be >= nwin"))
    haskey(NORM_FUNCS, norm) || throw(ArgumentError("Unknown spectrum_type: $norm."))

    win = _get_window(wintype[1], nwin, wintype[2] == :periodic ? true : false)
    frames = _get_frames(x, nwin, noverlap)
    wframes = _get_wframes(frames, win)

    nwin < nfft && (wframes = [vcat(wframe, zeros(Float64, nfft - nwin)) for wframe in wframes])

    plan = plan_rfft(1:nfft)
    spec = NORM_FUNCS[norm](combinedims(collect(plan * wframe for wframe in wframes)), win)
    halve && @views spec[[1, end], :] .*= 0.5 

    freq = range(0, stop = sr/2, length = size(spec, 1))

    Stft(sr, size(x, 1), StftSetup(nfft, nwin, noverlap, wintype, norm), StftData(spec, freq, win, frames))
end

function Base.show(io::IO, stft::Stft)
    print(io, "Stft(")
    print(io, "spec=$(size(stft.data.spec)), ")
    print(io, "nfft=$(stft.setup.nfft), ")
    print(io, "nwin=$(stft.setup.nwin), ")
    print(io, "noverlap=$(stft.setup.noverlap), ")
    print(io, "wintype=$(stft.setup.wintype), ")
    print(io, "norm=:$(stft.setup.norm))")
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
