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
#                                stft struct                                   #
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
	stft_freq :: Vector{Float64}
	info      :: NamedTuple

	function Stft{F}(
		stft_spec :: Matrix{T},
		stft_freq :: Vector{Float64},
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
#                                  utilities                                   #
#------------------------------------------------------------------------------#
function get_onesided_fft_range(stft_size::Int64)
	if mod(stft_size, 2) == 0
		return collect(1:Int(stft_size / 2 + 1))   # EVEN
	else
		return collect(1:Int((stft_size + 1) / 2))  # ODD
	end
end # get_onesided_fft_range

#------------------------------------------------------------------------------#
#                                   get stft                                   #
#------------------------------------------------------------------------------#
function get_stft(
	frames          :: AudioFrames;
	stft_size       :: Int64=get_wsize(frames),
	frequency_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr÷2),
	spectrum_type   :: Symbol=:power, # :power, :magnitude
)::Stft
	frequency_range isa Tuple && (frequency_range = FreqRange(first(frequency_range), last(frequency_range)))
	sr = get_info(frames).sr
	win_size   = get_wsize(frames)
	overlap    = get_ovrlap(frames)
	n_channels = nchannels(frames)
	f          = reduce(hcat, get_frames(frames))

    # validate frequency range
    (0 ≤ frequency_range.low < frequency_range.hi ≤ sr / 2) ||
        throw(ArgumentError("Frequency range must be (0, sr÷2). " *
				"Got range: $(frequency_range.low), $(frequency_range.hi), sr÷2 = $(sr÷2)"))
    
    # validate overlap
    (0 < overlap < win_size) ||
        throw(ArgumentError("Overlap length must be < window length. " *
				"Got overlap = $overlap, window length = $win_size"))
    
    # validate stft_size
    stft_size < win_size &&
        throw(ArgumentError("stft_size must be ≥ window length. " *
				"Got stft_size = $stft_size, window length = $win_size"))
    
    # Validate mono audio
    n_channels == 1 ||
        throw(ArgumentError("STFT is designed to work on mono audiofiles only. " *
				"Got $n_channels channels"))

	# ensure frames is of length stft_size
	# if the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
	# this zero-padding in the time domain results in an interpolation in the frequency domain, 
	# which can provide a more detailed view of the spectral content of the signal.
	@inline @views f = win_size < stft_size ? 
		vcat(f, zeros(eltype(frames), stft_size - win_size, size(f, 2))) : 
		f[1:stft_size, :]

	# get fft
	stft_spec = fft(f, (1,))

	# post process
	# trim to desired range
	bin_low  = ceil(Int,  frequency_range.low * stft_size / sr + 1)
	bin_high = floor(Int, frequency_range.hi  * stft_size / sr + 1)
	bins     = collect(bin_low:bin_high)
	stft_spec = @views stft_spec[bins, :]

	# convert to half-sided power or magnitude spectrum
	spectrum_funcs = Dict(
		:power     => frames -> real.(frames .* conj.(frames)),
		:magnitude => frames -> abs.(frames),
	)
	# check if spectrum_type is valid
	@assert haskey(spectrum_funcs, spectrum_type) "Unknown spectrum_type: $spectrum_type."

	stft_spec = spectrum_funcs[spectrum_type](stft_spec)

	# trim borders
	# halve the first bin if it's the lowest bin
	bin_low  == 1 && (@views stft_spec[1, :] *= 0.5)
	# halve the last bin if it's the Nyquist bin and FFT length is even
	bin_high == fld(stft_size, 2) + 1 && iseven(stft_size) && (@views stft_spec[end, :] *= 0.5)

	# create frequency vector
	stft_freq = (sr / stft_size) * (bins .- 1)
	# shift final bin if fftLength is odd and the final range is full to fs/2.
	if stft_size % 2 != 0 && bin_high == floor(fftLength / 2 + 1)
		stft_freq[end] = sr * (stft_size - 1) / (2 * stft_size)
	end

	info = (;
		sr,
		stft_size,
		win_size,
		overlap,
		frequency_range,
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
