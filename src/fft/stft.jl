# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct StftSetup <: AbstractSetup
    sr       :: Int64
    nfft     :: Int64
    winsize  :: Int64
    overlap  :: Int64
    spectrum :: Base.Callable
    window   :: Vector{<:Real}
end

# ---------------------------------------------------------------------------- #
#                                 stft struct                                  #
# ---------------------------------------------------------------------------- #
"""
    Stft{T} <: AbstractSpectrogram

A concrete implementation of `AbstractSpectrogram` that stores Short-Time Fourier Transform 
(STFT) data and associated metadata.

# Type Parameters
- `F`: Source type that generated this STFT (e.g., `Frames`)
- `T`: Element type of the spectral data matrix (e.g., `Float64`, `ComplexF64`)

# Fields
- `spec::Matrix{T}`: The STFT spectral data matrix where each column represents 
  a time frame and each row represents a frequency bin
- `freq::Vector{Float64}`: Frequency vector corresponding to the spectral bins (Hz)
- `info::NamedTuple`: Metadata containing analysis parameters including:
  - `sr`: Sample rate
  - `nfft`: FFT size used for analysis
  - `winsize`: Window size
  - `overlap`: Overlap between windows
  - `freqrange`: Frequency range analyzed
  - `spectrum`: Type of spectrum (`:power` or `:magnitude`)

# Constructor
    Stft{F}(spec::Matrix{T}, freq::Vector{Float64}, info::NamedTuple) where {T}

# Examples
```julia
# Get STFT from audio frames
frames = get_frames(audiofile)
stft   = get_stft(frames; nfft=512, spectrum=:power)

# Access data
spec_matrix = get_spec(stft)     # Get spectral data
frequencies = get_freq(stft)     # Get frequency vector
metadata    = get_setup(stft)     # Get analysis parameters
```

# See also:
[`AbstractSpectrogram`](@ref), [`get_spec`](@ref), [`get_freq`](@ref), [`get_setup`](@ref)
"""
struct Stft{T} <: AbstractSpectrogram
	spec :: Matrix{T}
	freq :: StepRangeLen
	info :: StftSetup

	function Stft(
		spec :: Matrix{T},
		freq :: StepRangeLen,
		info :: StftSetup
	) where {T<:AudioData}
		new{T}(spec, freq, info)
	end
end

#------------------------------------------------------------------------------#
#                                   methods                                    #
#------------------------------------------------------------------------------#
Base.eltype(::Stft{T}) where T = T

"""
    get_data(m::Stft{T}) -> Matrix

Get the stft spectrogram data matrix, transposed to (nframes × nbands).
"""
@inline get_data(s::Stft{T}) where T = s.spec

"""
    get_freq(m::Stft{T}) -> Vector

Get the frequency range in Hz.
"""
@inline get_freq(s::Stft{T}) where T = s.freq

"""
    get_setup(m::Stft{T}) -> StftSetup

Get the configuration metadata for the stft spectrogram.
"""
@inline get_setup(s::Stft{T}) where T = s.info

"""
    get_sr(s::Stft{T}) where T -> Int64

Get the sample rate used in STFT computation.
"""
@inline get_sr(s::Stft{T}) where T = s.info.sr

"""
    get_nfft(s::Stft{T}) where T -> Int64

Get the FFT size used in STFT computation.

The FFT size determines frequency resolution and the number of frequency bins.
For one-sided spectra (real signals), the number of bins is nfft/2 + 1.
"""
@inline get_nfft(s::Stft{T}) where T = s.info.nfft

"""
    get_spectrum(s::Stft{T}) where T -> Function

Get the spectrum type function used in STFT computation.
"""
@inline get_spectrum(s::Stft{T}) where T = s.info.spectrum

@inline get_winsize(s::Stft{T}) where T = s.info.winsize

@inline get_overlap(s::Stft{T}) where T = s.info.overlap

"""
    get_window(s::Stft{T}) where T -> Vector{<:Real}

Get the time-domain window coefficients used in STFT computation.
"""
@inline get_window(s::Stft{T}) where T = s.info.window

# ---------------------------------------------------------------------------- #
#                                     show                                     #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, s::Stft{T}) where T
    nfreqs, nframes = size(get_data(s))
    sr = s.info.sr
    spec_type = string(s.info.spectrum)
    
    print(io, "Stft{$T}(")
    print(io, "$nframes frames × $nfreqs bins, ")
    print(io, "sr=$sr Hz, ")
    print(io, "spectrum=$spec_type)")
end

function Base.show(io::IO, ::MIME"text/plain", s::Stft{T}) where T
    nframes = size(get_data(s), 2)
    sr = s.info.sr
    nfft = s.info.nfft
    winsize = s.info.winsize
    overlap = s.info.overlap
    spectrum_type = s.info.spectrum

    hop_size = winsize - overlap
    
    println(io, "Stft{$T}")
    println(io, "    Sample rate:     $sr Hz")
    println(io, "    Frames:          $nframes")
    println(io, "    FFT size:        $nfft")
    println(io, "    Window size:     $winsize samples")
    println(io, "    Overlap:         $overlap samples")
    println(io, "    Hop size:        $hop_size samples")
    println(io, "    Spectrum type:   $spectrum_type")
end

#------------------------------------------------------------------------------#
#                           spectrum normalizations                            #
#------------------------------------------------------------------------------#
"""
    power(f::AbstractArray) -> AbstractArray

Compute power spectral density from complex FFT values.

Transforms complex frequency-domain values to real power values by computing
the squared magnitude: |X(f)|² = real(X(f) · conj(X(f))).
"""
power(f)     = @. real(f * conj(f))

"""
    magnitude(f::AbstractArray) -> AbstractArray

Compute magnitude spectrum from complex FFT values.

Transforms complex frequency-domain values to real magnitude values by computing
the absolute value: |X(f)|.
"""
magnitude(f) = @. abs(f)

#------------------------------------------------------------------------------#
#                                  utilities                                   #
#------------------------------------------------------------------------------#
function _get_onesided_stft_range(nfft::Int64)::UnitRange
	return iseven(nfft) ?
		(1:nfft >> 1 + 1) : # even
        (1:(nfft + 1) >> 1) # odd
end

#------------------------------------------------------------------------------#
#                                   get stft                                   #
#------------------------------------------------------------------------------#
"""
    Stft(frames::Frames; nfft::Int64, spectrum::Base.Callable) -> Stft

Compute Short-Time Fourier Transform (STFT) from pre-computed audio frames.

Applies FFT to each windowed frame to obtain a time-frequency representation.
Zero-padding is automatically applied if `nfft > winsize` to achieve frequency 
interpolation without changing spectral content.

# Arguments
- `frames::Frames`: Pre-computed windowed audio frames

# Keyword Arguments
- `nfft::Int64`: FFT size (default: `winsize`). Must be ≥ window size.
  Larger values provide finer frequency resolution through interpolation.
- `spectrum::Base.Callable`: Spectrum type (default: `power`)
  - `power`: Power spectral density |X(f)|²
  - `magnitude`: Magnitude spectrum |X(f)|

# Returns
- `Stft`: STFT object containing spectral data and metadata

# Examples
```julia
# Basic usage with default FFT size
frames = Frames(audio; winsize=512, winstep=256)
stft = Stft(frames)

# With zero-padding for finer frequency resolution
stft = Stft(frames; nfft=1024)  # Double frequency bins

# Magnitude spectrum instead of power
stft = Stft(frames; spectrum=magnitude)
```

# Throws
- `ArgumentError`: If overlap ≥ winsize or nfft < winsize

# Notes
- Zero-padding in time domain = frequency interpolation
- Does not change spectral content, only interpolates between bins
- One-sided spectrum is returned (real signals are symmetric in frequency)
- Frequency resolution: Δf = sr / nfft

# See also
[`Frames`](@ref), [`power`](@ref), [`magnitude`](@ref)
"""
function Stft(
	frames   :: Frames;
	nfft     :: Int64=get_winsize(frames),
	spectrum :: Base.Callable=power, # power, magnitude
)::Stft
	sr        = get_setup(frames).sr
	winsize   = get_winsize(frames)
	overlap   = get_overlap(frames)
	winframes = get_winframes(frames)

    # validate overlap
    (0 < overlap < winsize) ||
        throw(ArgumentError("Overlap length must be < window length. " *
				"Got overlap = $overlap, window length = $winsize"))
    
    # validate nfft
    nfft < winsize &&
        throw(ArgumentError("nfft must be ≥ window length. " *
				"Got nfft = $nfft, window length = $winsize"))

	# ensure frames is of length nfft
	# if the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
	# this zero-padding in the time domain results in an interpolation in the frequency domain, 
	# which can provide a more detailed view of the spectral content of the signal.
    @inline @views winframes = winsize < nfft ? 
		vcat(winframes, zeros(eltype(winframes), nfft - winsize, size(winframes, 2))) : 
		winframes[1:nfft, :]

	# fft -> one side -> spectrum normalization
	spec = @views fft(winframes, (1,))[_get_onesided_stft_range(nfft), :] |> spectrum

	# frequency vector
	T = eltype(spec)
	freq = (0:size(spec, 1)-1) .* (T(sr) / T(nfft))

	info = StftSetup(
		sr,
		nfft,
		winsize,
		overlap,
		spectrum,
		get_window(frames)
	)

	return Stft(spec, freq, info)
end

"""
    Stft(audio::AudioFile; win::Base.Callable, type::Base.Callable, periodic::Bool, kwargs...) -> Stft

Compute Short-Time Fourier Transform (STFT) directly from an audio file.

Convenience constructor that combines frame extraction and STFT computation in one call.
Automatically selects appropriate default window parameters based on sample rate.

# Arguments
- `audio::AudioFile`: Input audio file

# Keyword Arguments
- `win::Base.Callable`: Window configuration function (default: adaptive based on sample rate)
  - For sr ≤ 8000 Hz: winsize=256, winstep=128
  - For sr > 8000 Hz: winsize=512, winstep=256
- `type::Base.Callable`: Window type function (default: `hanning`)
  - Examples: `hanning`, `hamming`, `blackman`, `rectangular`
- `periodic::Bool`: Use periodic window normalization (default: `true`)
- Additional kwargs are passed to `Stft(frames; kwargs...)`
  - `nfft::Int64`: FFT size
  - `spectrum::Base.Callable`: Spectrum type (`power` or `magnitude`)

# Returns
- `Stft`: STFT object containing spectral data and metadata

# Examples
```julia
# Default settings (adaptive based on sample rate)
audio = AudioFile("speech.wav")
stft = Stft(audio)

# Custom window configuration
stft = Stft(audio; 
    win=movingwindow(winsize=1024, winstep=512),
    type=hamming
)

# Custom FFT size and magnitude spectrum
stft = Stft(audio; 
    nfft=2048,
    spectrum=magnitude
)
```

# See also
[`Stft(::Frames)`](@ref), [`Frames`](@ref), [`AudioFile`](@ref), [`movingwindow`](@ref)
"""
function Stft(
	audio    :: AudioFile;
	winsize  :: Int64=get_sr(audio)≤8000 ? 256 : 512,
	winstep  :: Int64=get_sr(audio)≤8000 ? 128 : 256,
	type     :: Base.Callable=hanning,
    periodic :: Bool=true,
	kwargs...
	)::Stft
	frames = Frames(audio; winsize, winstep, type, periodic)
	Stft(frames; kwargs...)
end

