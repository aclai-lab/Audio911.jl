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

---

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

"""
function MelSpec(
    stft     :: Stft;
    win_norm :: Bool=true,
    kwargs...
)
    # validate scale if provided
    scale = get(kwargs, :scale, htk)
    scale in (htk, slaney) || throw(ArgumentError(
        "MelSpec only supports `htk` or `slaney` scale, got `$(nameof(scale))`."
    ))

    fbank = auditory_fbank(get_sr(stft); sfreq=get_freq(stft), nfft=get_nfft(stft), kwargs...)
    MelSpec(stft, fbank; win_norm)
end

# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
"""
    BarkSpec{F,B,T} <: AbstractSpectrogram

Bark-scale spectrogram of an audio signal.

`BarkSpec` is a mel-spectrogram variant in which the filterbank frequency axis
is warped according to the **Bark scale**, a psychoacoustic scale based on the
critical bands of human hearing. It models the non-linear frequency resolution
of the auditory system more faithfully in the mid-to-high frequency range than
the standard mel scale.

Internally, `BarkSpec` is backed by the same [`MelSpec`](@ref) / [`FBank`](@ref)
infrastructure, with the scale fixed to `bark`. All accessor methods are
forwarded accordingly.

# Fields
- `mel::MelSpec{F,B,T}`: Underlying mel spectrogram computed with `scale=bark`.

# Constructors

    BarkSpec(stft::Stft, fbank::AbstractFBank; win_norm::Bool=false)

Low-level constructor. Wraps a pre-computed [`FBank`](@ref) (must have been
designed with `scale=bark`).

---

    BarkSpec(stft::Stft; win_norm::Bool=true, kwargs...)

Convenience constructor. Automatically designs a Bark-scale filterbank from the
STFT parameters and computes the spectrogram.

## Arguments
- `stft::Stft`: Input short-time Fourier transform.

## Keywords
- `win_norm::Bool=true`: Apply window normalization (see [`MelSpec`](@ref)).
- `kwargs...`: Additional keyword arguments forwarded to [`auditory_fbank`](@ref),
  except `scale` which is always set to `bark`:
  - `nbands::Int64=26`: Number of Bark-spaced filter bands.
  - `norm::Function=bandwidth`: Filterbank normalization (`bandwidth`, `area`, `none_norm`).
  - `domain::Symbol=:linear`: Filter design domain (`:linear` or `:warped`).
  - `freqrange::FreqRange=(0, sr÷2)`: Frequency range in Hz.

# References
- Traunmüller, H. (1990). Analytical expressions for the tonotopic sensory scale.
  *Journal of the Acoustical Society of America*, 88(1), 97–100.
- Zwicker, E., & Fastl, H. (1999). *Psychoacoustics: Facts and Models* (2nd ed.).
  Springer.

# See also
[`MelSpec`](@ref), [`ErbSpec`](@ref), [`auditory_fbank`](@ref)
"""
struct BarkSpec{F,B,T} <: AbstractSpectrogram
    mel :: MelSpec{F,B,T}
end

function BarkSpec(
    stft     :: Stft,
    fbank    :: AbstractFBank;
    win_norm :: Bool=false
)
    BarkSpec(MelSpec(stft, fbank; win_norm))
end

function BarkSpec(
    stft     :: Stft;
    win_norm :: Bool=true,
    kwargs...
)
    # scale is fixed to bark, reject explicit override
    haskey(kwargs, :scale) && throw(ArgumentError(
        "BarkSpec does not accept a `scale` keyword; the Bark scale is always used. " *
        "Use MelSpec for htk/slaney scales."
    ))

    fbank = auditory_fbank(get_sr(stft); sfreq=get_freq(stft), nfft=get_nfft(stft), scale=bark, kwargs...)
    BarkSpec(stft, fbank; win_norm)
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    Base.eltype(::BarkSpec{F,B,T}) -> Type

Return the numeric element type `T` of the Bark spectrogram data.
"""
Base.eltype(::BarkSpec{F,B,T}) where {F,B,T} = T

"""
    get_data(b::BarkSpec) -> Matrix

Return the Bark spectrogram data matrix `(nframes × nbands)`.
"""
@inline get_data(b::BarkSpec)      = get_data(b.mel)

"""
    get_freq(b::BarkSpec) -> Vector

Return the center frequencies of the Bark filterbank bands, in Hz.
"""
@inline get_freq(b::BarkSpec)      = get_freq(b.mel)

"""
    get_freqrange(b::BarkSpec) -> FreqRange

Return the frequency range covered by the Bark filterbank, in Hz.
"""
@inline get_freqrange(b::BarkSpec) = get_freqrange(b.mel)

"""
    get_setup(b::BarkSpec) -> MelSpecSetup

Return the configuration metadata of the underlying [`MelSpec`](@ref).
"""
@inline get_setup(b::BarkSpec)     = get_setup(b.mel)

"""
    get_sr(b::BarkSpec) -> Int64

Return the sample rate (in Hz) of the source audio signal.
"""
@inline get_sr(b::BarkSpec)        = get_sr(b.mel)

"""
    get_nbands(b::BarkSpec) -> Int64

Return the number of Bark filterbank bands.
"""
@inline get_nbands(b::BarkSpec)    = get_nbands(b.mel)

# ---------------------------------------------------------------------------- #
#                                     show                                     #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, b::BarkSpec{F,B,T}) where {F,B,T}
    nframes, nbands = size(get_data(b))
    print(io, "BarkSpec{$F,$B,$T}($nframes frames × $nbands bands, sr=$(get_sr(b)) Hz)")
end

function Base.show(io::IO, ::MIME"text/plain", b::BarkSpec{F,B,T}) where {F,B,T}
    nframes, nbands = size(get_data(b))
    freq_range = extrema(get_freq(b))
    println(io, "BarkSpec{$F,$B,$T}")
    println(io, "  Dimensions:")
    println(io, "    Frames:     $nframes")
    println(io, "    Bark bands: $nbands")
    println(io, "  Configuration:")
    println(io, "    Sample rate:       $(get_sr(b)) Hz")
    println(io, "    Frequency range:   $(round(freq_range[1], digits=1)) - $(round(freq_range[2], digits=1)) Hz")
    print(io,   "    Window normalized: $(get_setup(b).win_norm)")
end

# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
"""
    ErbSpecSetup <: AbstractSetup

Configuration structure for ERB (gammatone) spectrogram computation.

# Fields
- `sr::Int64`: Sample rate of the source audio signal, in Hz.
- `win_norm::Bool`: Whether window normalization was applied.

# See also
[`ErbSpec`](@ref), [`gammatone_fbank`](@ref)
"""
struct ErbSpecSetup <: AbstractSetup
    sr       :: Int64
    win_norm :: Bool
end

# ---------------------------------------------------------------------------- #
#                               erb spectrogram                                #
# ---------------------------------------------------------------------------- #
"""
    ErbSpec{F,B,T} <: AbstractSpectrogram

ERB-scale spectrogram of an audio signal using a gammatone filterbank.

`ErbSpec` models the frequency selectivity of the human cochlea by applying a
bank of **gammatone filters** spaced uniformly on the **Equivalent Rectangular
Bandwidth (ERB)** scale. Compared to triangular mel or Bark filters, gammatone
filters provide a more physiologically accurate model of auditory processing:

- Smooth, asymmetric bandpass shapes matching basilar membrane responses.
- ERB spacing closely follows the frequency resolution of the inner ear.
- Well-suited for auditory modelling, hearing aid research, and biologically
  inspired speech/audio analysis.

`ErbSpec` does **not** accept a `scale` keyword; the ERB scale and gammatone
filter shape are intrinsic to this representation. For triangular-filter
alternatives use [`MelSpec`](@ref) (htk/slaney) or [`BarkSpec`](@ref) (bark).

# Fields
- `spec::Matrix{T}`: Spectrogram matrix `(nbands × nframes)`.
- `fbank::FBank`: Gammatone filterbank used for computation.
- `info::ErbSpecSetup`: Configuration metadata (sample rate, window norm flag).

# Constructors

    ErbSpec(stft::Stft, fbank::AbstractFBank; win_norm::Bool=false)

Low-level constructor. Applies a pre-computed gammatone [`FBank`](@ref) to `stft`.

---

    ErbSpec(stft::Stft; win_norm::Bool=true, kwargs...)

Convenience constructor. Automatically designs an ERB-scale gammatone filterbank
from the STFT parameters and computes the spectrogram.

## Arguments
- `stft::Stft`: Input short-time Fourier transform.

## Keywords
- `win_norm::Bool=true`: Apply window normalization to account for the analysis
  window's energy (see [`MelSpec`](@ref) for details).
- `kwargs...`: Additional keyword arguments forwarded to [`gammatone_fbank`](@ref):
  - `nbands::Int64=26`: Number of ERB-spaced gammatone filters.
  - `norm::Function=bandwidth`: Filter normalization (`bandwidth`, `area`, `none_norm`).
  - `freqrange::FreqRange=(0, sr÷2)`: Frequency range in Hz.

# References
- Patterson, R. D., Nimmo-Smith, I., Holdsworth, J., & Rice, P. (1988).
  An efficient auditory filterbank based on the gammatone function.
  *APU Report 2341*, MRC Applied Psychology Unit, Cambridge.
- Slaney, M. (1993). *An Efficient Implementation of the Patterson-Holdsworth
  Auditory Filter Bank*. Apple Computer Technical Report 35.
- Ellis, D. P. W. (2009). Gammatone-like spectrograms.
  [https://labrosa.ee.columbia.edu/matlab/gammatonegram](https://labrosa.ee.columbia.edu/matlab/gammatonegram)

# See also
[`MelSpec`](@ref), [`BarkSpec`](@ref), [`gammatone_fbank`](@ref)
"""
struct ErbSpec{F,B,T} <: AbstractSpectrogram
    spec  :: Matrix{T}
    fbank :: FBank
    info  :: ErbSpecSetup

    function ErbSpec{F,B}(
        spec  :: Matrix{T},
        fbank :: FBank,
        info  :: ErbSpecSetup
    ) where {F,T,B}
        new{F,B,T}(spec, fbank, info)
    end
end

function ErbSpec(
    stft     :: Stft,
    fbank    :: AbstractFBank;
    win_norm :: Bool=false
)::ErbSpec
    spectrum   = get_spectrum(stft)
    window     = get_window(stft)
    spec       = get_data(stft)
    filterbank = get_data(fbank)

    if win_norm
        filterbank = copy(filterbank)
        spectrum == power     && (filterbank *= winpower(1, window))
        spectrum == magnitude && (filterbank *= winmagnitude(1, window))
    end

    info = ErbSpecSetup(get_sr(stft), win_norm)

    ErbSpec{typeof(stft),typeof(fbank)}(filterbank * spec, fbank, info)
end

function ErbSpec(
    stft     :: Stft;
    win_norm :: Bool=true,
    kwargs...
)
    fbank = gammatone_fbank(get_sr(stft); nfft=get_nfft(stft), kwargs...)
    ErbSpec(stft, fbank; win_norm)
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    Base.eltype(::ErbSpec{F,B,T}) -> Type

Return the numeric element type `T` of the ERB spectrogram data.
"""
Base.eltype(::ErbSpec{F,B,T}) where {F,B,T} = T

"""
    get_data(e::ErbSpec) -> Matrix

Return the ERB spectrogram data matrix `(nframes × nbands)`.
"""
@inline get_data(e::ErbSpec)      = e.spec'

"""
    get_freq(e::ErbSpec) -> Vector

Return the center frequencies of the gammatone filterbank bands, in Hz.
"""
@inline get_freq(e::ErbSpec)      = get_freq(e.fbank)

"""
    get_freqrange(e::ErbSpec) -> FreqRange

Return the frequency range covered by the gammatone filterbank, in Hz.
"""
@inline get_freqrange(e::ErbSpec) = get_freqrange(e.fbank)

"""
    get_setup(e::ErbSpec) -> ErbSpecSetup

Return the [`ErbSpecSetup`](@ref) configuration metadata.
"""
@inline get_setup(e::ErbSpec)     = e.info

"""
    get_sr(e::ErbSpec) -> Int64

Return the sample rate (in Hz) of the source audio signal.
"""
@inline get_sr(e::ErbSpec)        = e.info.sr

"""
    get_nbands(e::ErbSpec) -> Int64

Return the number of gammatone filterbank bands.
"""
@inline get_nbands(e::ErbSpec)    = get_nbands(e.fbank)

# ---------------------------------------------------------------------------- #
#                                     show                                     #
# ---------------------------------------------------------------------------- #
function Base.show(io::IO, e::ErbSpec{F,B,T}) where {F,B,T}
    nframes, nbands = size(get_data(e))
    print(io, "ErbSpec{$F,$B,$T}($nframes frames × $nbands bands, sr=$(get_sr(e)) Hz)")
end

function Base.show(io::IO, ::MIME"text/plain", e::ErbSpec{F,B,T}) where {F,B,T}
    nframes, nbands = size(get_data(e))
    freq_range = extrema(get_freq(e))
    println(io, "ErbSpec{$F,$B,$T}")
    println(io, "  Dimensions:")
    println(io, "    Frames:     $nframes")
    println(io, "    ERB bands:  $nbands")
    println(io, "  Configuration:")
    println(io, "    Sample rate:       $(get_sr(e)) Hz")
    println(io, "    Frequency range:   $(round(freq_range[1], digits=1)) - $(round(freq_range[2], digits=1)) Hz")
    print(io,   "    Window normalized: $(e.info.win_norm)")
end
