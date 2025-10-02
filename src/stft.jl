# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
"""
    AbstractSpectrogram

Abstract supertype for all spectrogram representations in the Audio911.jl package.

This type serves as the base for concrete spectrogram implementations that store
time-frequency representations of audio signals. Subtypes should implement
methods for accessing spectral data, frequency information, and metadata.

# Concrete Subtypes
- [`Stft`](@ref): Short-Time Fourier Transform spectrogram representation

# Interface
Concrete subtypes are expected to provide:
- Spectral data matrix (typically complex or real-valued)
- Frequency vector corresponding to spectral bins
- Metadata about the analysis parameters

# See also
[`Stft`](@ref), [`get_stft`](@ref), [`get_freq`](@ref), [`get_info`](@ref)
"""
abstract type AbstractSpectrogram end

#------------------------------------------------------------------------------#
#                            spectrogram methods                               #
#------------------------------------------------------------------------------#
"""
    get_spec(s::AbstractSpectrogram) -> Matrix{T}

Extract the spectrogram matrix from an `AbstractSpectrogram` container.

# See also: [`get_freq`](@ref), [`get_info`](@ref)
"""
get_spec(s::AbstractSpectrogram) = error("get_spec is not implemented for type $(typeof(s)).")

"""
    get_freq(s::AbstractSpectrogram) -> Vector{Float64}

Extract the frequency vector from an `AbstractSpectrogram` container.

# See also: [`get_spec`](@ref), [`get_info`](@ref)
"""
get_freq(s::AbstractSpectrogram) = error("get_freq is not implemented for type $(typeof(s)).")

"""
    get_info(s::AbstractSpectrogram) -> NamedTuple

Extract metadata from an `AbstractSpectrogram` container.

# See also: [`get_spec`](@ref), [`get_freq`](@ref)
"""
get_info(s::AbstractSpectrogram) = error("get_info is not implemented for type $(typeof(s)).")

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
- `stft_spec::Matrix{T}`: The STFT spectral data matrix where each column represents 
  a time frame and each row represents a frequency bin
- `stft_freq::Vector{Float64}`: Frequency vector corresponding to the spectral bins (Hz)
- `info::NamedTuple`: Metadata containing analysis parameters including:
  - `sr`: Sample rate
  - `stft_size`: FFT size used for analysis
  - `win_size`: Window size
  - `overlap`: Overlap between windows
  - `frequency_range`: Frequency range analyzed
  - `spectrum_type`: Type of spectrum (`:power` or `:magnitude`)

# Constructor
    Stft{F}(stft_spec::Matrix{T}, stft_freq::Vector{Float64}, info::NamedTuple) where {F,T}

# Examples
```julia
# Get STFT from audio frames
frames = get_frames(audiofile)
stft = get_stft(frames; stft_size=512, spectrum_type=:power)

# Access data
spec_matrix = get_spec(stft)     # Get spectral data
frequencies = get_freq(stft)     # Get frequency vector
metadata    = get_info(stft)     # Get analysis parameters
```

# See also:
[`AbstractSpectrogram`](@ref), [`get_spec`](@ref), [`get_freq`](@ref), [`get_info`](@ref)
"""
struct Stft{F,T} <: AbstractSpectrogram
	stft_spec :: Matrix{T}
	stft_freq :: StepRangeLen
	info      :: NamedTuple

	function Stft{F}(
		stft_spec :: Matrix{T},
		stft_freq :: StepRangeLen,
		info      :: NamedTuple
	) where {F,T}
		new{F,T}(stft_spec, stft_freq, info)
	end
end

#------------------------------------------------------------------------------#
#                                 stft methods                                 #
#------------------------------------------------------------------------------#
Base.eltype(::Stft{T}) where T = T

get_spec(s::Stft) = s.stft_spec
get_freq(s::Stft) = s.stft_freq
get_info(s::Stft) = s.info

#------------------------------------------------------------------------------#
#                           spectrum normalizations                            #
#------------------------------------------------------------------------------#
power(f) = @. real(f * conj(f))
magnitude(f) = @. abs(f)

#------------------------------------------------------------------------------#
#                                  utilities                                   #
#------------------------------------------------------------------------------#
function get_onesided_fft_range(stft_size::Int64)::UnitRange
	return iseven(stft_size) ?
		(1:stft_size >> 1 + 1) : # even
        (1:(stft_size + 1) >> 1) # odd
end

#------------------------------------------------------------------------------#
#                                   get stft                                   #
#------------------------------------------------------------------------------#
function get_stft(
	frames          :: AudioFrames;
	stft_size       :: Int64=get_winsize(frames),
	# frequency_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr÷2),
	spectrum_type   :: Base.Callable=power, # power, magnitude
)::Stft
	# frequency_range isa Tuple && (frequency_range = FreqRange(first(frequency_range), last(frequency_range)))
	sr = get_info(frames).sr
	win_size   = get_winsize(frames)
	overlap    = get_overlap(frames)
	winframes  = get_winframes(frames)

    # validate overlap
    (0 < overlap < win_size) ||
        throw(ArgumentError("Overlap length must be < window length. " *
				"Got overlap = $overlap, window length = $win_size"))
    
    # validate stft_size
    stft_size < win_size &&
        throw(ArgumentError("stft_size must be ≥ window length. " *
				"Got stft_size = $stft_size, window length = $win_size"))

	# ensure frames is of length stft_size
	# if the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
	# this zero-padding in the time domain results in an interpolation in the frequency domain, 
	# which can provide a more detailed view of the spectral content of the signal.
    @inline @views winframes = win_size < stft_size ? 
		vcat(winframes, zeros(eltype(winframes), stft_size - win_size, size(winframes, 2))) : 
		winframes[1:stft_size, :]

	# fft -> one side -> spectrum normalization
	stft_spec = @views fft(winframes, (1,))[get_onesided_fft_range(stft_size), :] |> spectrum_type

	# frequency vector
	stft_freq = (0:size(stft_spec, 1)-1) .* (sr / stft_size)

	info = (;
		sr,
		stft_size,
		win_size,
		overlap,
		spectrum_type,
	)
	info = merge(info, get_info(frames))

	return Stft{AudioFrames}(stft_spec, stft_freq, info)
end

function get_stft(
	afile :: AudioFile;
    win   :: WinFunction=MovingWindow(
                            window_size=samplerate(afile)≤8000 ? 256 : 512,
                            window_step=samplerate(afile)≤8000 ? 128 : 256
                        ),
	type  :: Tuple{Symbol, Symbol}=(:hann, :periodic),
	kwargs...
	)::Stft
	frames = get_frames(afile; win, type)
	get_stft(frames; kwargs...)
end
