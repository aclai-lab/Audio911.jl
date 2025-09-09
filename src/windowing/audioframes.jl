# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
"""
    AbstractWinFunction

Base type for windowing function implementations.
"""
abstract type AbstractWinFunction end

"""
    AbstractAudioFrames

Abstract base type for audio frame containers.
"""
abstract type AbstractAudioFrames end

# ---------------------------------------------------------------------------- #
#                                win functions                                 #
# ---------------------------------------------------------------------------- #
"""
    WinFunction <: AbstractWinFunction

Callable wrapper for windowing algorithms with parameters.

# Fields
- `func::Function`: The windowing implementation function
- `params::NamedTuple`: Algorithm-specific parameters
"""
struct WinFunction <: AbstractWinFunction
    func   :: Function
    params :: NamedTuple
end
# Make it callable - npoints is passed at execution time
(w::WinFunction)(npoints::Int; kwargs...) = w.func(npoints; w.params..., kwargs...)

"""
    MovingWindow(; window_size::Int, window_step::Int) -> WinFunction

Create a moving window that slides across the time series.

# Parameters
- `window_size`: Number of time points in each window
- `window_step`: Step size between consecutive windows

# Example
    win = MovingWindow(window_size=10, window_step=5)
    intervals = win(100)  # For 100-point time series
"""
function MovingWindow(;
    window_size::Int,
    window_step::Int,
)
    WinFunction(movingwindow, (;window_size, window_step))
end

# ---------------------------------------------------------------------------- #
#                                 AudioFrames                                  #
# ---------------------------------------------------------------------------- #
"""
    AudioFrames{T} <: AbstractAudioFrames

Container for windowed audio data with associated windowing parameters.

This struct holds the result of applying a windowing function to audio data,
storing both the windowed frames and the parameters used to generate them.

# Type Parameters
- `T`: The numeric type of the audio data elements (e.g., `Float32`, `Float64`)

# Fields
- `frames::Vector{<:AudioFormat{T}}`: Vector of windowed audio frames, where each frame 
  is either a `Vector{T}` (mono) or `Matrix{T}` (multi-channel)
- `win::WinFunction`: The windowing function used to generate the frames 
  (e.g., `MovingWindow`, `AdaptiveWindow`, `SplitWindow`, `WholeWindow`)
- `type::Tuple{Symbol, Symbol}`: Window type specification as `(window_shape, periodicity)`,
  where `window_shape` can be `:hann`, `:hamming`, `:blackman`, etc., and `periodicity` 
  can be `:periodic` or `:symmetric`

# Constructor
    AudioFrames(frames, win, type)

# Example
```julia
# Load audio file
audiofile = load("audio.wav")

# Create windowing function
win = AdaptiveWindow(nwindows=3, relative_overlap=0.1)

# Generate frames
frames = get_frames(audiofile; win=win, type=(:hann, :periodic))

# Access the windowed data
windowed_data = frames.frames
window_function = frames.win
window_type = frames.type
```

# See Also
- [`get_frames`](@ref): Function to create `AudioFrames` from audio data
- [`WinFunction`](@ref): Base type for windowing functions
- [`AudioFormat`](@ref): Type alias for audio data formats
"""
struct AudioFrames{T} <: AbstractAudioFrames
    frames :: Vector{<:AudioFormat{T}}
    win    :: WinFunction
	type   :: Tuple{Symbol, Symbol}
    info   :: NamedTuple

    function AudioFrames(
        frames :: Vector{<:AudioFormat{T}},
        win    :: WinFunction,
        type   :: Tuple{Symbol, Symbol},
        info   :: NamedTuple
    ) where T
    new{T}(frames, win, type, info)
    end
end

#------------------------------------------------------------------------------#
#                                    methods                                   #
#------------------------------------------------------------------------------#
Base.length(f::AudioFrames) = length(f.frames)
Base.eltype(::AudioFrames{T}) where T = T

"""
    get_frames(f::AudioFrames) -> Vector{<:AudioFormat}

Extract the windowed audio frames from an `AudioFrames` container.

# See Also
- [`get_frames(::AudioFile)`](@ref): Function to create `AudioFrames` from audio data
- [`AudioFrames`](@ref): Container type for windowed audio data
"""
get_frames(f::AudioFrames) = f.frames

"""
    nchannels(f::AudioFrames) -> Int

Get the number of audio channels from an `AudioFrames` container.

# See Also
- [`AudioFrames`](@ref): Container type for windowed audio data
- [`get_frames`](@ref): Function to create audio frames
- [`Base.length`](@ref): Get number of frames in the container
"""
nchannels(f::AudioFrames) = size(f.frames[1], 2)

"""
    get_wsize(f::AudioFrames) -> Int

Get the window size parameter from an `AudioFrames` object.

# See Also
- [`get_wstep`](@ref): Get window step size
- [`get_ovrlap`](@ref): Get overlap length
"""
get_wsize(f::AudioFrames)  = f.win.params.window_size

"""
    get_wstep(f::AudioFrames) -> Int

Get the window step parameter from an `AudioFrames` object.

# See Also
- [`get_wsize`](@ref): Get window size
- [`get_ovrlap`](@ref): Get overlap length
"""
get_wstep(f::AudioFrames)  = f.win.params.window_step

"""
    get_ovrlap(f::AudioFrames) -> Int

Get the overlap length from an `AudioFrames` object.

# See Also
- [`get_wsize`](@ref): Get window size  
- [`get_wstep`](@ref): Get window step size
"""
get_ovrlap(f::AudioFrames) = f.win.params.window_size - f.win.params.window_step

"""
    get_info(f::AudioFrames) -> NamedTuple

Extract metadata from an `AudioFrames` container.

# See Also
- [`AudioFrames`](@ref): Container type for windowed audio data
- [`get_frames`](@ref): Function to create audio frames
- [`get_wsize`](@ref), [`get_wstep`](@ref), [`get_ovrlap`](@ref): Access specific windowing parameters
"""
get_info(f::AudioFrames)   = f.info

#------------------------------------------------------------------------------#
#                                    frames                                    #
#------------------------------------------------------------------------------#
"""
    get_frames(X; win=AdaptiveWindow(nwindows=3, relative_overlap=0.1), type=(:hann, :periodic)) -> AudioFrames

Apply windowing to an audio file and return windowed frames with applied window functions.

This function segments the audio data according to the specified windowing strategy,
applies the chosen window function to each frame, and returns an `AudioFrames` container
with the processed data and metadata.

# Arguments
- `X::AudioFile`: Input audio file containing the audio data to be windowed

# Keyword Arguments
- `win::WinFunction`: Windowing strategy to use. Default is `AdaptiveWindow(nwindows=3, relative_overlap=0.1)`
  - `AdaptiveWindow`: Creates overlapping windows with adaptive sizing
  - `MovingWindow`: Sliding window with fixed size and step
  - `SplitWindow`: Equal non-overlapping segments
  - `WholeWindow`: Single window encompassing entire audio
- `type::Tuple{Symbol, Symbol}`: Window function specification as `(shape, periodicity)`. 
  Default is `(:hann, :periodic)`
  - First element (shape): `:hann`, `:hamming`, `:blackman`, `:bartlett`, etc.
  - Second element (periodicity): `:periodic` or `:symmetric`

# Returns
- `AudioFrames{T}`: Container holding the windowed frames, windowing function, and parameters

# Examples
```julia
# Load audio file
audiofile = load("speech.wav")

# Basic usage with default parameters
frames = get_frames(audiofile)

# Using moving window with Hamming window
frames = get_frames(audiofile; 
                   win=MovingWindow(window_size=512, window_step=256),
                   type=(:hamming, :symmetric))

# Split into 4 equal segments with Blackman window
frames = get_frames(audiofile;
                   win=SplitWindow(nwindows=4),
                   type=(:blackman, :periodic))

# Adaptive windowing with high overlap
frames = get_frames(audiofile;
                   win=AdaptiveWindow(nwindows=5, relative_overlap=0.5),
                   type=(:hann, :periodic))

# Access the results
windowed_data = frames.frames      # Vector of windowed audio segments
window_func = frames.win           # The windowing function used
window_type = frames.type          # Window shape and periodicity
```

# Notes
- Each frame in the output has the window function applied element-wise
- Frame lengths may vary depending on the windowing strategy used
- The window function is automatically sized to match each frame's length
- All frames are returned as vectors (mono audio) or matrices (multi-channel)

# See Also
- [`AudioFrames`](@ref): Container type for windowed audio data
- [`AdaptiveWindow`](@ref), [`MovingWindow`](@ref), [`SplitWindow`](@ref), [`WholeWindow`](@ref): Windowing strategies
- [`WinFunction`](@ref): Base windowing function type
"""
function get_frames(
	afile :: AudioFile;
    win   :: WinFunction=MovingWindow(
                            window_size=sr(afile) ≤ 8000 ? 256 : 512,
                            window_step=sr(afile) ≤ 8000 ? 128 : 256
                        ),
	type  :: Tuple{Symbol, Symbol}=(:hann, :periodic),
)::AudioFrames
    # run the windowing algo and set windows indexes
    intervals = win(length(afile))

    frames = map(intervals) do interval
        frame      = nchannels(afile) == 1 ? data(afile)[interval] : data(afile)[interval, :]
        window, _  = gencoswin(type[1], length(interval), type[2])
        frame .* window
    end

    info = (;
        win_func = win.func,
        win_size = win.params.window_size,
        win_step = win.params.window_step,
        win_type = type
    )
    return AudioFrames(frames, win, type, info)
end

# function logEnergyCoeffs(
#         x::AbstractArray{T}
# ) where {T <: Real}
#     DT = eltype(x)
#     E = sum(x .^ 2, dims = 1) # eleva tutti gli elementi ^2 e li somma per colonna
#     E[E .== 0] .= floatmin(DT) # se un valore è zero, lo sostituisce col valore più piccolo positivo possibile, in accordo col tipo utilizzato
#     logE = log.(E) # fa il log di tutti gli elementi
# end # logEnergyCoeffs

# function windowing(
#         x::Union{AbstractVector{T}, AbstractArray{T}},
#         fftLength::Int64 = 256,
#         winType::Symbol = :hann,
#         winParam::Symbol = :symmetric,
#         logEnergy::Bool = false
# ) where {T <: Real}
#     xLength = size(x, 1) # lunghezza audio
#     nChan = size(x, 2) # numero canali (mono, stereo)
#     DT = eltype(x) # restituisce il tipo degli elementi

#     # parto con un if then ma sarebbe bello implementare un Dict
#     if (winType == :hann || winType == :hamming || winType == :blackman ||
#         winType == :flattopwin)
#         win, hL = gencoswin(winType, fftLength, winParam)
#     elseif (winType == :rect)
#         win, hL = rectwin(fftLength)
#     end

#     wincast = convert.(DT, win[:]) # casta al tipo originale del file audio
#     nHops = Integer(floor((xLength - fftLength) / hL) + 1) # numero delle window necessarie
#     y = buffer(x, fftLength, hL)

#     # calcola il vettore coefficienti energia. in matlab lo fa qui, ma forse è logico spostarlo dopo il windowing?
#     if (logEnergy == true)
#         logE = logEnergyCoeffs(y)
#     end

#     return y' .* wincast'
# end # function windowing

# function fade(
#         x::Union{AbstractVector{T}, AbstractArray{T}},
#         fftLength::Int64,
#         type::Symbol
# ) where {T <: Real}
#     xLength = size(x, 1) # lunghezza audio
#     nChan = size(x, 2) # numero canali (mono, stereo)
#     DT = eltype(x) # restituisce il tipo degli elementi

#     win, hL = gencoswin(:hann, fftLength, :symmetric)

#     wincast = convert.(DT, win[:]) # casta al tipo originale del file audio

#     # for c = 1:numChan
#     if (type == :in)
#         for w in 1:Int(round(fftLength / 2))
#             x[w] *= wincast[w]
#         end
#     elseif (type == :out)
#         for w in Int(round(fftLength / 2)):-1:fftLength
#             x[end - fftLength + w] *= wincast[w]
#         end
#     end
#     # end

#     return x
# end # function fade

#------------------------------------------------------------------------------#
#                                   windowing                                  #
#------------------------------------------------------------------------------#
# function _get_frames(
# 	x::AbstractVector{Float64};
# 	window_type::Tuple{Symbol, Symbol},
# 	window_length::Int64,
# 	overlap_length::Int64,
# )
# 	frames = buffer(x, window_length, window_length - overlap_length)
# 	window, _ = gencoswin(window_type[1], window_length, window_type[2])

# 	return frames, window
# end

# function _get_frames(x::AbstractVector{Float64}, s::AudioSetup)
# 	_get_frames(x, window_type = s.window_type, window_length = s.window_length, overlap_length = s.overlap_length)
# end

# function get_frames(
# 	x::AbstractVector{<:AbstractFloat},
# 	window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# 	window_length::Int64 = 256,
# 	overlap_length::Int64 = 128,
# )
# 	@assert 0 < overlap_length < window_length "Overlap length must be < window length."

# 	frames, window = _get_frames(
# 		eltype(x) == Float64 ? x : Float64.(x),
# 		window_type = window_type,
# 		window_length = window_length,
# 		overlap_length = overlap_length,
# 	)
# end

# get_frames!(a::AudioObj; kwargs...) = a.data.frames, a.setup.window = get_frames(a.data.x; kwargs...)

# function _get_frames2(x::AbstractVector{Float64}, s::AudioSetup)
# 	_get_frames2(x, window_type = s.window_type, window_length = s.window_length, overlap_length = s.overlap_length)
# end

# function get_frames2(
# 	x::AbstractVector{<:AbstractFloat},
# 	window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# 	window_length::Int64 = 256,
# 	overlap_length::Int64 = 128,
# )
# 	@assert 0 < overlap_length < window_length "Overlap length must be < window length."

# 	frames = _get_frames2(
# 		eltype(x) == Float64 ? x : Float64.(x),
# 		window_type = window_type,
# 		window_length = window_length,
# 		overlap_length = overlap_length,
# 	)
# end

# get_frames2!(a::AudioObj; kwargs...) = a.data.frames, a.setup.window = get_frames2(a.data.x; kwargs...)
