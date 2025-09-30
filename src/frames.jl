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
#                                    types                                     #
# ---------------------------------------------------------------------------- #
# availables window functions from package DSP
const AVAIL_WINDOWS = (
    rect, hanning, hamming, cosine, lanczos, triang,
    bartlett, bartlett_hann, blackman
)

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
    get_size(w::WinFunction) -> Int

Get the window size parameter from a `WinFunction` object.

# See also: [`get_step`](@ref), [`WinFunction`](@ref), [`MovingWindow`](@ref)
"""
get_size(w::WinFunction) = w.params.window_size

"""
    get_step(w::WinFunction) -> Int

Get the window step parameter from a `WinFunction` object.

# See also: [`get_size`](@ref), [`WinFunction`](@ref), [`MovingWindow`](@ref): Moving window constructor
"""
get_step(w::WinFunction) = w.params.window_step

"""
    get_overlap(w::WinFunction) -> Int

Get the overlap length between consecutive windows from a `WinFunction` object.

# See also: [`get_size`](@ref), [`get_step`](@ref), [`WinFunction`](@ref), [`MovingWindow`](@ref)
"""
get_overlap(w::WinFunction) = get_size(w) - get_step(w)

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
)::WinFunction
    WinFunction(movingwindow, (;window_size, window_step))
end

# ---------------------------------------------------------------------------- #
#                                 AudioFrames                                  #
# ---------------------------------------------------------------------------- #
"""
    AudioFrames{T} <: AbstractAudioFrames

Container for windowed audio data with associated windowing parameters.

This struct holds the result of applying a windowing function to audio data,
storing the windowed frames as a matrix where each column represents one frame,
along with the window function and metadata used to generate them.

# Type Parameters
- `T`: The numeric type of the audio data elements (e.g., `Float32`, `Float64`)

# Fields
- `frames::Matrix{T}`: Matrix of windowed audio frames where each column is one frame.
  For mono audio, each column contains the windowed samples for that frame.
  For multi-channel audio, the matrix structure preserves channel information.
- `window::Vector{Float64}`: The actual window function values that were applied to each frame
- `info::NamedTuple`: Metadata containing:
  - `sr::Int`: Sample rate of the original audio
  - `win::WinFunction`: The windowing function object used
  - `type::Symbol`: window type function from DSP package

# Constructor
    AudioFrames(frames::Vector{<:AudioFormat{T}}, window::Vector{Float64}, info::NamedTuple)

The constructor automatically converts a vector of individual frames into a matrix format
where each column represents one frame, using `reduce(hcat, frames)` for efficient storage.

# Matrix Structure
- **Rows**: Sample indices within each frame (window_size samples per frame)
- **Columns**: Frame indices (number of frames depends on audio length and windowing parameters)

# Example
```julia
# Load audio file
audiofile = load("audio.wav")

# Create windowing function
win = MovingWindow(window_size=512, window_step=256)

# Generate frames with Hann window
frames = get_frames(audiofile; win=win, type=hanning)

# Access the data
frames_matrix = frames.frames        # Matrix where each column is a frame
window_values = frames.window        # Window function values applied
sample_rate = frames.info.sr         # Original sample rate
window_func = frames.info.win        # Windowing function used
```

# See also: [`get_frames`](@ref), [`WinFunction`](@ref), [`MovingWindow`](@ref), [`AudioFormat`](@ref)
"""
struct AudioFrames{T} <: AbstractAudioFrames
    frames :: Matrix{T}
    window :: Vector{T}
    info   :: NamedTuple

    function AudioFrames(
        frames :: Vector{<:AudioFormat{T}},
        window :: Vector{Float64},
        info   :: NamedTuple
    ) where T
        frames_matrix = reduce(hcat, frames)
        new{T}(frames_matrix, window, info)
    end

    function AudioFrames(
        frames :: Matrix{<:AudioFormat{T}},
        window :: Vector{Float64},
        info   :: NamedTuple
    ) where T
        new{T}(frames, window, info)
    end
end

#------------------------------------------------------------------------------#
#                                    methods                                   #
#------------------------------------------------------------------------------#
Base.length(f::AudioFrames) = size(f.frames, 2)
Base.eltype(::AudioFrames{T}) where T = T

"""
    get_frames(f::AudioFrames) -> Vector{<:AudioFormat}

Extract the buffered audio frames from an `AudioFrames` container.

# See also: [`get_frames(::AudioFile)`](@ref), [`AudioFrames`](@ref)
"""
get_frames(f::AudioFrames) = f.frames

"""
    get_window(f::AudioFrames) -> Vector{Float64}

Extract the window from an `AudioFrames` container.

# See also: [`get_frames`](@ref), [`get_wframes`](@ref), [`AudioFrames`](@ref)
"""
get_window(f::AudioFrames) = f.window

"""
    get_wframes(f::AudioFrames) -> Matrix

Get the windowed frames with the window function applied element-wise.

# See also: [`get_frames`](@ref), [`get_window`](@ref), [`AudioFrames`](@ref)
"""
get_winframes(f::AudioFrames) = get_frames(f) .* get_window(f)

"""
    nchannels(f::AudioFrames) -> Int

Get the number of audio channels from an `AudioFrames` container.

# See also: [`AudioFrames`](@ref), [`get_frames`](@ref), [`Base.length`](@ref)
"""
nchannels(f::AudioFrames) = size(f.frames[1], 2)

"""
    get_size(f::AudioFrames) -> Int

Get the window size parameter from an `AudioFrames` object.

# See also: [`get_wstep`](@ref), [`get_overlap`](@ref)
"""
get_size(f::AudioFrames)  = f.info.win.params.window_size

"""
    get_wstep(f::AudioFrames) -> Int

Get the window step parameter from an `AudioFrames` object.

# See also: [`get_size`](@ref), [`get_overlap`](@ref)
"""
get_step(f::AudioFrames)  = f.info.win.params.window_step

"""
    get_overlap(f::AudioFrames) -> Int

Get the overlap length from an `AudioFrames` object.

# See also: [`get_size`](@ref), [`get_step`](@ref)
"""
get_overlap(f::AudioFrames) = get_size(f) - get_step(f)

"""
    get_info(f::AudioFrames) -> NamedTuple

Extract metadata from an `AudioFrames` container.

# See also: [`AudioFrames`](@ref)
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
- `type::Base.Callable`: 

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

# See also: [`AudioFrames`](@ref), [`WinFunction`](@ref), [`MovingWindow`](@ref)
"""
function get_frames(
	afile    :: AudioFile;
    win      :: WinFunction=MovingWindow(
                window_size=samplerate(afile) ≤ 8000 ? 256 : 512,
                window_step=samplerate(afile) ≤ 8000 ? 128 : 256
            ),
	type     :: Base.Callable=hanning,
    periodic :: Bool=true
)::AudioFrames
    # setup parameters
    afile_length = length(afile)

    # check invalid parameters
    afile_length < win.params.window_size && throw(ArgumentError("Audio file length ($afile_length)" * 
        " is shorter than window size ($(win.params.window_size))"))
    in(type, AVAIL_WINDOWS) || throw(ArgumentError("Window type $(type) not supported." * 
        " Available windows: $(AVAIL_WINDOWS)"))

    # buffer audio file
    intervals = win(afile_length)
    frames = map(intervals) do interval
        nchannels(afile) == 1 ? data(afile)[interval] : data(afile)[interval, :]
    end
    @show size(frames,1)
    @show size(frames,2)
    @show length(frames[1])

    # get window from DSP package
    win_size = get_size(win)
    window = if periodic
        # zerophase and circshift are used for compatibility with Matlab's `periodic` windows.
        circshift(type(win_size; zerophase=true), -(win_size >> 1))
    else
        type(win_size)
    end

    # collect infos
    info = (;
        sr   = samplerate(afile),
        win  = win,
        type = type
    )

    return AudioFrames(frames, window, info)
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
