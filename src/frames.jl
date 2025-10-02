# ---------------------------------------------------------------------------- #
#                               abstract type                                  #
# ---------------------------------------------------------------------------- #
"""
    AbstractFrame

Base type for windowing function implementations.
"""
abstract type AbstractFrame end

# ---------------------------------------------------------------------------- #
#                              abstract methods                                #
# ---------------------------------------------------------------------------- #
"""
    get_size(f::AbstractFrame) -> Int

Get the window size parameter from a `AbstractFrame` object.

# See also: [`AbstractFrame`](@ref), [`get_step`](@ref), [`get_overlap`](@ref)
"""
get_size(f::AbstractFrame) = error("get_size is not implemented for type $(typeof(f)).")

"""
    get_step(f::AbstractFrame) -> Int

Get the window step parameter from a `AbstractFrame` object.

# See also: [`AbstractFrame`](@ref), [`get_size`](@ref), [`get_overlap`](@ref)
"""
get_step(f::AbstractFrame) = error("get_step is not implemented for type $(typeof(f)).")

"""
    get_overlap(f::AbstractFrame) -> Int

Get the overlap length between consecutive windows from a `AbstractFrame` object.

# See also: [`AbstractFrame`](@ref), [`get_size`](@ref), [`get_step`](@ref)
"""
get_overlap(f::AbstractFrame) = error("get_overlap is not implemented for type $(typeof(f)).")

"""
    get_frames(f::AbstractFrame) -> Vector{<:AudioFormat}

Extract the buffered audio frames from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref), [`get_winframes`](@ref), [`get_window`](@ref)
"""
get_frames(f::AbstractFrame) = error("get_frames is not implemented for type $(typeof(f)).")

"""
    get_window(f::AbstractFrame) -> Vector{Float64}

Extract the window from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref), [`get_frames`](@ref), [`get_winframes`](@ref)
"""
get_window(f::AbstractFrame) = error("get_window is not implemented for type $(typeof(f)).")

"""
    get_winframes(f::AbstractFrame) -> Matrix

Get the windowed frames with the window function applied element-wise.

# See also: [`AbstractFrame`](@ref), [`get_frames`](@ref), [`get_window`](@ref)
"""
get_winframes(f::AbstractFrame) = error("get_winframes is not implemented for type $(typeof(f)).")

"""
    get_winsize(f::AbstractFrame) -> Int

Get the length of the window function from an `AbstractFrame` container.

# See also: [`AudioFrames`](@ref), [`get_window`](@ref), [`get_winframes`](@ref)
"""
get_winsize(f::AbstractFrame) = error("get_winsize is not implemented for type $(typeof(f)).")

"""
    get_info(f::AbstractFrame) -> NamedTuple

Extract metadata from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref)
"""
get_info(f::AbstractFrame) = error("get_info is not implemented for type $(typeof(f)).")

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
    WinFunction <: AbstractFrame

Callable wrapper for windowing algorithms with parameters.

# Fields
- `func::Function`    : The windowing implementation function
- `params::NamedTuple`: Algorithm-specific parameters
"""
struct WinFunction <: AbstractFrame
    func   :: Function
    params :: NamedTuple
end

# Make it callable - npoints is passed at execution time
(w::WinFunction)(npoints::Int; kwargs...) = w.func(npoints; w.params..., kwargs...)

# ---------------------------------------------------------------------------- #
#                                   methods                                    #
# ---------------------------------------------------------------------------- #
get_size(w::WinFunction)    = w.params.window_size
get_step(w::WinFunction)    = w.params.window_step
get_overlap(w::WinFunction) = get_size(w) - get_step(w)

# ---------------------------------------------------------------------------- #
#                                 constructor                                  #
# ---------------------------------------------------------------------------- #
"""
    MovingWindow(; window_size::Int, window_step::Int) -> WinFunction

Create a moving window function that slides across the time series.

# Parameters
- `window_size`: Number of time points in each window
- `window_step`: Step size between consecutive windows

# Example
    win = MovingWindow(window_size=10, window_step=5)
    intervals = win(100)  # For 100-point time series
"""
function MovingWindow(;
    size::Int,
    step::Int,
)::WinFunction
    WinFunction(movingwindow, (window_size=size, window_step=step))
end

# ---------------------------------------------------------------------------- #
#                                 AudioFrames                                  #
# ---------------------------------------------------------------------------- #
"""
    AudioFrames{T} <: AbstractFrame

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
struct AudioFrames{T} <: AbstractFrame
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
function AudioFrames(
    audiofile :: AudioFormat,
    sr        :: Int64;
    win       :: WinFunction=MovingWindow(
        size=sr ≤ 8000 ? 256 : 512,
        step=sr ≤ 8000 ? 128 : 256
    ),
	type      :: Base.Callable=hanning,
    periodic  :: Bool=true
)
    # auto mono conversion
    size(audiofile, 2) > 1 && (audiofile = convert2mono(audiofile))
    
    # setup parameters
    file_length = length(audiofile)

    # check invalid parameters
    file_length < get_size(win) && throw(ArgumentError("Audio file length ($file_length)" * 
        " is shorter than window size ($(get_size(win)))"))
    in(type, AVAIL_WINDOWS) || throw(ArgumentError("Window type $(type) not supported." * 
        " Available windows: $(AVAIL_WINDOWS)"))

    # buffer audio file
    intervals = win(file_length)
    frames = [audiofile[i] for i in intervals]

    # get window from DSP package
    win_size = get_size(win)
    window = if periodic
        # zerophase and circshift are used for compatibility with Matlab's `periodic` windows.
        circshift(type(win_size; zerophase=true), -(win_size >> 1))
    else
        type(win_size)
    end

    # collect infos
    info = (; sr, win, type)

    return AudioFrames(frames, window, info)
end

AudioFrames(a::AudioFile; kwargs...) = AudioFrames(data(a), samplerate(a); kwargs...)

#------------------------------------------------------------------------------#
#                                    methods                                   #
#------------------------------------------------------------------------------#
Base.length(f::AudioFrames) = size(f.frames, 2)
Base.eltype(::AudioFrames{T}) where T = T

get_size(f::AudioFrames)  = f.info.win.params.window_size
get_step(f::AudioFrames)  = f.info.win.params.window_step
get_overlap(f::AudioFrames) = get_size(f) - get_step(f)

get_frames(f::AudioFrames) = f.frames
get_window(f::AudioFrames) = f.window
get_winframes(f::AudioFrames) = get_frames(f) .* get_window(f)
get_info(f::AudioFrames)   = f.info
