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
    Stft{F,T} <: AbstractSpectrogram

A concrete implementation of `AbstractSpectrogram` that stores Short-Time Fourier Transform 
(STFT) data and associated metadata.

# Type Parameters
- `F`: Source type that generated this STFT (e.g., `AudioFrames`)
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
  - `freq_range`: Frequency range analyzed
  - `spectrum`: Type of spectrum (`:power` or `:magnitude`)

# Constructor
    Stft{F}(spec::Matrix{T}, freq::Vector{Float64}, info::NamedTuple) where {F,T}

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
struct Stft{F,T} <: AbstractSpectrogram
	spec :: Matrix{T}
	freq :: StepRangeLen
	info :: StftSetup

	function Stft{F}(
		spec :: Matrix{T},
		freq :: StepRangeLen,
		info :: StftSetup
	) where {F,T}
		new{F,T}(spec, freq, info)
	end
end

#------------------------------------------------------------------------------#
#                           spectrum normalizations                            #
#------------------------------------------------------------------------------#
power(f)     = @. real(f * conj(f))
magnitude(f) = @. abs(f)

#------------------------------------------------------------------------------#
#                                  utilities                                   #
#------------------------------------------------------------------------------#
function get_onesided_stft_range(nfft::Int64)::UnitRange
	return iseven(nfft) ?
		(1:nfft >> 1 + 1) : # even
        (1:(nfft + 1) >> 1) # odd
end

#------------------------------------------------------------------------------#
#                                   get stft                                   #
#------------------------------------------------------------------------------#
function Stft(
	frames   :: AudioFrames;
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
	spec = @views fft(winframes, (1,))[get_onesided_stft_range(nfft), :] |> spectrum

	# frequency vector
	freq = (0:size(spec, 1)-1) .* (sr / nfft)

	info = StftSetup(
		sr,
		nfft,
		winsize,
		overlap,
		spectrum,
		get_window(frames)
	)

	return Stft{AudioFrames}(spec, freq, info)
end

function Stft(
	afile    :: AudioFile;
    win      :: Base.Callable=movingwindow(
		winsize=get_sr(afile)≤8000 ? 256 : 512,
		winstep=get_sr(afile)≤8000 ? 128 : 256
	),
	type     :: Base.Callable=hanning,
    periodic :: Bool=true,
	kwargs...
	)::Stft
	frames = AudioFrames(afile; win, type, periodic)
	Stft(frames; kwargs...)
end

#------------------------------------------------------------------------------#
#                                 stft methods                                 #
#------------------------------------------------------------------------------#
Base.eltype(::Stft{T}) where T = T

@inline get_data(s::Stft) = s.spec
@inline get_freq(s::Stft) = s.freq
@inline get_setup(s::Stft) = s.info

@inline get_sr(s::Stft)       = s.info.sr
@inline get_nfft(s::Stft)     = s.info.nfft
@inline get_spectype(s::Stft) = s.info.spectrum
@inline get_window(s::Stft)   = s.info.window
