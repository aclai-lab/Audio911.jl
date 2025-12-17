# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct MelSpecSetup <: AbstractSetup
	sr       :: Int64
	win_norm :: Bool
end

# ---------------------------------------------------------------------------- #
#                          linear spectrogram struct                           #
# ---------------------------------------------------------------------------- #
struct MelSpec{F,B,T} <: AbstractSpectrogram
    spec  :: Matrix{T}
	fbank :: FBank
    info  :: MelSpecSetup

    function MelSpec{F,B}(
        spec  :: Matrix{T},
		fbank :: FBank,
        info  :: MelSpecSetup
    ) where {F,T,B}
        new{F,B,T}(spec, fbank, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    Base.eltype(::MelSpec{T}) -> Type

Return the element type of the mel spectrogram data.
"""
Base.eltype(::MelSpec{F,B,T}) where {F,B,T} = T

"""
    get_data(m::MelSpec) -> Matrix

Get the mel spectrogram data matrix, transposed to (nframes × nbands).
"""
@inline get_data(m::MelSpec)  = m.spec'

"""
    get_freq(m::MelSpec) -> Vector

Get the frequency bins in Hz.
"""
@inline get_freq(m::MelSpec)  = get_freq(m.fbank)

"""
    get_freqrange(m::MelSpec) -> Vector

Get the frequency range in Hz.
"""
@inline get_freqrange(m::MelSpec)  = get_freqrange(m.fbank)

"""
    get_setup(m::MelSpec) -> MelSpecSetup

Get the configuration metadata for the mel spectrogram.
"""
@inline get_setup(m::MelSpec) = m.info

"""
    get_sr(s::MelSpec) -> Int64

Get the sample rate used in mel spectrogram computation.
"""
@inline get_sr(m::MelSpec) = m.info.sr

"""
    get_nbands(s::MelSpec) -> Int64

Get the number of filter bands in the filterbank.
"""
@inline get_nbands(m::MelSpec) = get_nbands(m.fbank)

# ---------------------------------------------------------------------------- #
#                                     show                                     #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, m::MelSpec{F,B,T}) where {F,B,T}
    nframes, nbands = size(get_data(m))
    sr = m.info.sr
    win_norm = m.info.win_norm ? "enabled" : "disabled"
    
    print(io, "MelSpec{$F,$B,$T}(")
    print(io, "$nframes frames × $nbands bands, ")
    print(io, "sr=$sr Hz, ")
    print(io, "win_norm=$win_norm)")
end

function Base.show(io::IO, ::MIME"text/plain", m::MelSpec{F,B,T}) where {F,B,T}
    nframes, nbands = size(get_data(m))
    sr = m.info.sr
    win_norm = m.info.win_norm
    
    freq_range = extrema(get_freq(m))
    duration = nframes / sr  # approximate duration
    
    println(io, "MelSpec{$F,$B,$T}")
    println(io, "  Dimensions:")
    println(io, "    Frames:     $nframes")
    println(io, "    Mel bands:  $nbands")
    println(io, "  Configuration:")
    println(io, "    Sample rate:        $sr Hz")
    println(io, "    Frequency range:    $(round(freq_range[1], digits=1)) - $(round(freq_range[2], digits=1)) Hz")
    println(io, "    Window normalized:  $win_norm")
    print(io,   "    Duration (approx):  $(round(duration, digits=3)) s")
end

# ---------------------------------------------------------------------------- #
#                                 get mel spec                                 #
# ---------------------------------------------------------------------------- #
"""
    MelSpec(stft::Stft, fbank::AbstractFBank; win_norm::Bool=false) -> MelSpec

Compute mel spectrogram from STFT using a provided filterbank.

Applies the mel filterbank to the STFT to transform from linear frequency bins 
to mel-frequency bands. Optionally applies window normalization to account for 
the analysis window's energy.

# Arguments
- `stft::Stft`: Short-time Fourier transform to process
- `fbank::AbstractFBank`: Mel filterbank for frequency transformation

# Keywords
- `win_norm::Bool=false`: Apply window normalization factor

# Examples
```julia
audiofile = Audio911.load(wav_file, format=Float64)
# Create STFT from audio frames
frames = Frames(audiofile; win=movingwindow(winsize=512, winstep=256))
stft = Stft(frames; spectrum=power)

# Create mel filterbank
fbank = auditory_fbank(16000; sfreq=get_freq(stft), nbands=40)

# Compute mel spectrogram
mel = MelSpec(stft, fbank; win_norm=true)
```
"""
function MelSpec(
	stft     :: Stft,
    fbank    :: AbstractFBank;
	win_norm :: Bool=false
)::MelSpec
    spectrum   = get_spectrum(stft)
	window     = get_window(stft)
	spec       = get_data(stft)
	filterbank = get_data(fbank)

	if win_norm
		spectrum == power     && (filterbank *= winpower(1, window))
		spectrum == magnitude && (filterbank *= winmagnitude(1, window))
	end

	info = MelSpecSetup(get_sr(stft), win_norm)

	return MelSpec{typeof(stft),typeof(fbank)}(filterbank * spec, fbank, info)
end

"""
    MelSpec(stft::Stft; win_norm::Bool=true, kwargs...) -> MelSpec

Compute mel spectrogram from STFT with automatic filterbank creation.

Convenience constructor that automatically creates an appropriate mel filterbank 
based on the STFT parameters before computing the mel spectrogram.

# Arguments
- `stft::Stft`: Short-time Fourier transform to process

# Keywords
- `win_norm::Bool=true`: Apply window normalization
- `kwargs...`: Additional arguments passed to `auditory_fbank`:
  - `freqrange::Tuple`: Frequency range in Hz (default: (0, sr/2))
  - `nbands::Int`: Number of mel bands (default: 40)
  - `norm::Symbol`: Filterbank normalization (`:bandwidth`, `:area`, `:none`)
  - `domain::Symbol`: Design domain (`:linear`, `:mel`)
  - `scale::Symbol`: Mel scale type (`:htk`, `:slaney`)

# Returns
- `MelSpec`: Mel spectrogram

# Examples
```julia
frames = AudioFrames(audio; win=movingwindow(winsize=512, winstep=256))
stft = Stft(frames; spectrum=power)

# With default mel filterbank
mel = MelSpec(stft)

# With custom filterbank parameters
mel = MelSpec(stft; 
    win_norm=true, 
    freqrange=(100, 8000),
    nbands=26,
    norm=:bandwidth,
    scale=:htk
)
```
"""
function MelSpec(
	stft     :: Stft;
    win_norm :: Bool=true,
    kwargs... # auditory filterbank kwargs
)
    fbank = auditory_fbank(get_sr(stft); sfreq=get_freq(stft), nfft=get_nfft(stft), kwargs...)
    MelSpec(stft, fbank; win_norm)
end
