# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct LinSpecSetup <: AbstractSetup
    sr            :: Int64
    freqrange     :: FreqRange
    spectrum_type :: Base.Callable
	win_norm      :: Bool
end

# ---------------------------------------------------------------------------- #
#                          linear spectrogram struct                           #
# ---------------------------------------------------------------------------- #
struct LinSpec{F,T} <: AbstractSpectrogram
    spec :: Matrix{T}
    freq :: StepRangeLen
    info :: LinSpecSetup

    function LinSpec{F}(
        spec :: Matrix{T},
        freq :: StepRangeLen,
        info :: LinSpecSetup
    ) where {F,T}
        new{F,T}(spec, freq, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
Base.eltype(::LinSpec{T}) where T = T

"""
    get_data(m::LinSpec) -> Matrix

Get the linear spectrogram data matrix, transposed to (nframes × nbands).
"""
get_data(s::LinSpec)  = s.spec'

"""
    get_freq(m::LinSpec) -> Vector

Get the frequency range in Hz.
"""
get_freq(s::LinSpec)  = s.freq

"""
    get_setup(m::LinSpec) -> MelSpecSetup

Get the configuration metadata for the linear spectrogram.
"""
get_setup(s::LinSpec) = s.info

# ---------------------------------------------------------------------------- #
#                                     show                                     #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, s::LinSpec{F,T}) where {F,T}
    nfreqs, nframes = size(get_data(s))
    freq_range = extrema(get_freq(s))
    sr = s.info.sr
    win_norm = s.info.win_norm ? "enabled" : "disabled"
    
    print(io, "LinSpec{$F,$T}(")
    print(io, "$nframes frames × $nfreqs bins, ")
    print(io, "$(round(freq_range[1], digits=0))-$(round(freq_range[2], digits=0)) Hz, ")
    print(io, "sr=$sr Hz, ")
    print(io, "win_norm=$win_norm)")
end

function Base.show(io::IO, ::MIME"text/plain", s::LinSpec{F,T}) where {F,T}
    nfreqs, nframes = size(get_data(s))
    freq = get_freq(s)
    info = s.info
    
    freq_range = extrema(freq)
    freq_step = step(freq)
    duration = nframes / info.sr  # approximate duration
    
    println(io, "LinSpec{$F,$T}")
    println(io, "  Dimensions:")
    println(io, "    Frames:          $nframes")
    println(io, "    Frequency bins:  $nfreqs")
    println(io, "  Configuration:")
    println(io, "    Sample rate:        $(info.sr) Hz")
    println(io, "    Frequency range:    $(round(freq_range[1], digits=1)) - $(round(freq_range[2], digits=1)) Hz")
    println(io, "    Frequency step:     $(round(freq_step, digits=2)) Hz")
    println(io, "    Spectrum type:      $(info.spectrum_type)")
    println(io, "    Window normalized:  $(info.win_norm)")
    print(io,   "    Duration (approx):  $(round(duration, digits=3)) s")
end

# ---------------------------------------------------------------------------- #
#                                  utilities                                   #
# ---------------------------------------------------------------------------- #
function get_freq_range(
    freqrange :: FreqRange,
    stft_size  :: Int64,
    sr         :: Int64
)
    # convert frequencies to bin indices
    bin_low  = cld(get_low(freqrange) * stft_size, sr) + 1
    bin_high = fld(get_hi(freqrange)  * stft_size, sr) + 1

    return bin_low, bin_high
end

# ---------------------------------------------------------------------------- #
#                                 get lin spec                                 #
# ---------------------------------------------------------------------------- #
"""
    LinSpec(stft::Stft; freqrange::FreqRange=(0, sr÷2), win_norm::Bool=false) -> LinSpec

Compute linear frequency spectrogram from Short-Time Fourier Transform (STFT).

Creates a linear spectrogram by extracting magnitude or power spectrum from STFT frames,
optionally restricted to a specific frequency range and with window normalization applied.

# Arguments
- `stft::Stft`: Short-Time Fourier Transform to process

# Keywords
- `freqrange::FreqRange=(0, sr÷2)`: Frequency range to extract in Hz as tuple `(low, high)`
  - Default: Full Nyquist range from 0 Hz to half the sample rate
  - Allows focusing on specific frequency bands (e.g., `(100, 8000)` for speech)
- `win_norm::Bool=false`: Apply window normalization to compensate for analysis window energy
  - When `true`, applies normalization factor based on spectrum type and window function
  - Accounts for power/magnitude loss due to windowing in time domain

# Examples

```julia
# Load audio and create frames
audio = Audio911.load("speech.wav", format=Float64)
frames = Frames(audio; win=movingwindow(winsize=512, winstep=256))

# Compute STFT with power spectrum
stft = Stft(frames; spectrum=power)

# Default: full frequency range, no normalization
linspec = LinSpec(stft)

# Focus on speech frequency range (100-8000 Hz)
linspec_speech = LinSpec(stft; freqrange=(100, 8000))

# Apply window normalization for accurate power measurements
linspec_norm = LinSpec(stft; win_norm=true)

# Access data
spec_data = get_data(linspec)  # (nframes × nfreqs) matrix
freq_axis = get_freq(linspec)  # frequency vector in Hz
```
"""
function LinSpec(
	stft      :: Stft;
	freqrange :: FreqRange=(0, get_setup(frames).sr>>1),
	win_norm  :: Bool=false
)::LinSpec
	spec = get_data(stft)
	freq = get_freq(stft)

	sr            = get_sr(stft)
	stft_size     = get_nfft(stft)
	spectrum_type = get_spectrum(stft)
	window        = get_window(stft)

	if freqrange != (0, sr >> 1)
		bin_low, bin_high = get_freq_range(freqrange, stft_size, sr)
		spec = @views spec[bin_low:bin_high, :]
		freq = freq[bin_low:bin_high]
	end

	win_norm_func = eval(Symbol("win" * string(spectrum_type)))
	win_norm && (spec = win_norm_func(spec, window))

	info = LinSpecSetup(sr, freqrange, spectrum_type, win_norm)

	return LinSpec{typeof(stft)}(spec .* 2, freq, info)
end

