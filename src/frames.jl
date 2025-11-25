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
    get_data(f::AbstractFrame) -> Vector{<:AudioFormat}

Extract the buffered audio frames from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref), [`get_winframes`](@ref), [`get_window`](@ref)
"""
get_data(f::AbstractFrame) = error("get_data is not implemented for type $(typeof(f)).")

"""
    get_window(f::AbstractFrame) -> Vector{Float64}

Extract the window from an `AbstractFrame` container.

# See also: [`AbstractFrame`](@ref), [`get_data`](@ref), [`get_winframes`](@ref)
"""
get_window(f::AbstractFrame) = error("get_window is not implemented for type $(typeof(f)).")

"""
    get_winframes(f::AbstractFrame) -> Matrix

Get the windowed frames with the window function applied element-wise.

# See also: [`AbstractFrame`](@ref), [`get_data`](@ref), [`get_window`](@ref)
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
#                                     info                                     #
# ---------------------------------------------------------------------------- #
struct FramesInfo <: AbstractInfo
    sr   :: Int64
    win  :: Base.Callable
    type :: Base.Callable
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
frames = get_data(audiofile; win=win, type=hanning)

# Access the data
frames_matrix = frames.frames        # Matrix where each column is a frame
window_values = frames.window        # Window function values applied
sample_rate = frames.info.sr         # Original sample rate
window_func = frames.info.win        # Windowing function used
```

# See also: [`get_data`](@ref), [`WinFunction`](@ref), [`MovingWindow`](@ref), [`AudioFormat`](@ref)
"""
struct AudioFrames{T} <: AbstractFrame
    frames :: Matrix{T}
    window :: Vector{T}
    info   :: FramesInfo

    function AudioFrames(
        frames :: Vector{<:AudioFormat{T}},
        window :: Vector{Float64},
        info   :: FramesInfo
    ) where T
        frames_matrix = reduce(hcat, frames)
        new{T}(frames_matrix, window, info)
    end

    function AudioFrames(
        frames :: Matrix{<:AudioFormat{T}},
        window :: Vector{Float64},
        info   :: FramesInfo
    ) where T
        new{T}(frames, window, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    frames                                    #
# ---------------------------------------------------------------------------- #
"""
    get_data(X; win=AdaptiveWindow(nwindows=3, relative_overlap=0.1), type=(:hann, :periodic)) -> AudioFrames

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
frames = get_data(audiofile)

# Using moving window with Hamming window
frames = get_data(audiofile; 
                   win=MovingWindow(window_size=512, window_step=256),
                   type=(:hamming, :symmetric))

# Split into 4 equal segments with Blackman window
frames = get_data(audiofile;
                   win=SplitWindow(nwindows=4),
                   type=(:blackman, :periodic))

# Adaptive windowing with high overlap
frames = get_data(audiofile;
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

# See also: [`AudioFrames`](@ref)
"""
function AudioFrames(
    audiofile :: AudioFormat,
    sr        :: Int64;
    win       :: Base.Callable=movingwindow(
        winsize=sr ≤ 8000 ? 256 : 512,
        winstep=sr ≤ 8000 ? 128 : 256
    ),
	type      :: Base.Callable=hanning,
    periodic  :: Bool=true
)
    # auto mono conversion
    size(audiofile, 2) > 1 && (audiofile = _convert_mono(audiofile))
    
    # setup parameters
    file_length = length(audiofile)

    # check invalid parameters
    file_length < win.winsize && throw(ArgumentError("Audio file length ($file_length)" * 
        " is shorter than window size ($(get_size(win)))"))
    in(type, AVAIL_WINDOWS) || throw(ArgumentError("Window type $(type) not supported." * 
        " Available windows: $(AVAIL_WINDOWS)"))

    # buffer audio file
    intervals = win(file_length)
    frames = [audiofile[i] for i in intervals]

    # get window from DSP package
    window = if periodic
        # zerophase and circshift are used for compatibility with Matlab's `periodic` windows.
        circshift(type(win.winsize; zerophase=true), -(win.winsize >> 1))
    else
        type(win.winsize)
    end

    # collect infos
    info = FramesInfo(sr, win, type)

    return AudioFrames(frames, window, info)
end

AudioFrames(a::AudioFile; kwargs...) = AudioFrames(get_data(a), get_sr(a); kwargs...)

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
Base.length(f::AudioFrames) = size(f.frames, 2)
Base.eltype(::AudioFrames{T}) where T = T

get_size(f::AudioFrames)    = f.info.win.winsize
get_step(f::AudioFrames)    = f.info.win.winstep
get_overlap(f::AudioFrames) = get_size(f) - get_step(f)

get_data(f::AudioFrames)      = f.frames
get_window(f::AudioFrames)    = f.window
get_winsize(f::AudioFrames)   = length(f.window)
get_winframes(f::AudioFrames) = get_data(f) .* get_window(f)
get_info(f::AudioFrames)      = f.info

# ---------------------------------------------------------------------------- #
#                                   base.show                                  #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, ::MIME"text/plain", f::AudioFrames{T}) where T
    n_frames = length(f)
    frame_size = get_size(f)
    step_size = get_step(f)
    overlap = get_overlap(f)
    sr = f.info.sr
    win_type = f.info.type
    
    # Calculate total duration
    total_samples = frame_size + (n_frames - 1) * step_size
    duration = total_samples / sr
    
    println(io, "AudioFrames{$T}")
    println(io, "  Frames:      $n_frames")
    println(io, "  Frame size:  $frame_size samples")
    println(io, "  Step:        $step_size samples")
    println(io, "  Overlap:     $overlap samples ($(round(100 * overlap / frame_size, digits=1))%)")
    println(io, "  Window:      $win_type")
    println(io, "  Sample rate: $sr Hz")
end

function Base.show(io::IO, f::AudioFrames{T}) where T
    n_frames = length(f)
    frame_size = get_size(f)
    print(io, "AudioFrames{$T}($n_frames frames × $frame_size samples)")
end
