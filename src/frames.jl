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
struct FramesSetup <: AbstractSetup
    sr   :: Int64
    win  :: Base.Callable
    type :: Base.Callable
end

# ---------------------------------------------------------------------------- #
#                                    Frames                                    #
# ---------------------------------------------------------------------------- #
struct Frames{T} <: AbstractFrame
    frames :: Matrix{T}
    window :: Vector{T}
    info   :: FramesSetup

    function Frames(
        frames :: Vector{<:AudioFormat{T}},
        window :: Vector{Float64},
        info   :: FramesSetup
    ) where T
        frames_matrix = reduce(hcat, frames)
        new{T}(frames_matrix, window, info)
    end

    function Frames(
        frames :: Matrix{<:AudioFormat{T}},
        window :: Vector{Float64},
        info   :: FramesSetup
    ) where T
        new{T}(frames, window, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
Base.length(f::Frames{T}) where T = size(f.frames, 2)
Base.eltype(::Frames{T})  where T = T

get_size(f::Frames{T}) where T = f.info.win.winsize
get_step(f::Frames{T}) where T = f.info.win.winstep
get_overlap(f::Frames{T}) where T = get_size(f) - get_step(f)

get_data(f::Frames{T}) where T = f.frames
get_window(f::Frames{T}) where T = f.window
get_winsize(f::Frames{T}) where T = length(f.window)
get_winframes(f::Frames{T}) where T = get_data(f) .* get_window(f)
get_setup(f::Frames{T}) where T = f.info
get_sr(f::Frames{T}) where T = f.info.sr

# ---------------------------------------------------------------------------- #
#                                   base.show                                  #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, ::MIME"text/plain", f::Frames{T}) where T
    n_frames   = length(f)
    frame_size = get_size(f)
    step_size  = get_step(f)
    overlap    = get_overlap(f)
    sr         = f.info.sr
    win_type   = f.info.type
    
    println(io, "Frames{$T}")
    println(io, "  Sample rate: $sr Hz")
    println(io, "  Frames:      $n_frames")
    println(io, "  Frame size:  $frame_size samples")
    println(io, "  Step:        $step_size samples")
    println(io, "  Overlap:     $overlap samples ($(round(100 * overlap / frame_size, digits=1))%)")
    println(io, "  Window:      $win_type")
end

function Base.show(io::IO, f::Frames{T}) where T
    n_frames   = length(f)
    frame_size = get_size(f)
    print(io, "Frames{$T}($n_frames frames × $frame_size samples)")
end

# ---------------------------------------------------------------------------- #
#                                    frames                                    #
# ---------------------------------------------------------------------------- #
"""
    Frames(X; win=adaptivewindow(nwindows=3, overlap=0.1), type=(:hann, :periodic)) -> Frames

Apply windowing to an audio file and return windowed frames with applied window functions.

This function segments the audio data according to the specified windowing strategy,
applies the chosen window function to each frame, and returns an `Frames` container
with the processed data and metadata.

# Arguments
- `X::AudioFile`: Input audio file containing the audio data to be windowed

# Keyword Arguments
- `win::WinFunction`: Windowing strategy to use. Default is `adaptivewindow(nwindows=3, overlap=0.1)`
  - `adaptivewindow`: Creates overlapping windows with adaptive sizing
  - `MovingWindow`: Sliding window with fixed size and step
  - `SplitWindow`: Equal non-overlapping segments
  - `WholeWindow`: Single window encompassing entire audio
- `type::Base.Callable`: 

# Returns
- `Frames{T}`: Container holding the windowed frames, windowing function, and parameters

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
                   win=adaptivewindow(nwindows=5, overlap=0.5),
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
"""
function Frames(
    audiofile :: AudioFormat,
    sr        :: Int64;
    winsize   :: Int64=sr ≤ 8000 ? 256 : 512,
    winstep   :: Int64=sr ≤ 8000 ? 128 : 256,
	type      :: Base.Callable=hanning,
    periodic  :: Bool=true
)
    # auto mono conversion
    size(audiofile, 2) > 1 && (audiofile = _convert_mono(audiofile))
    
    # setup parameters
    file_length = length(audiofile)

    # check invalid parameters
    file_length < winsize && throw(ArgumentError("Audio file length ($file_length)" * 
        " is shorter than window size ($(get_size(win)))"))
    in(type, AVAIL_WINDOWS) || throw(ArgumentError("Window type $(type) not supported." * 
        " Available windows: $(AVAIL_WINDOWS)"))

    # buffer audio file
    win = DataTreatments.movingwindow(; winsize, winstep)
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
    info = FramesSetup(sr, win, type)

    return Frames(frames, window, info)
end

Frames(a::AudioFile; kwargs...) = Frames(get_data(a), get_sr(a); kwargs...)

