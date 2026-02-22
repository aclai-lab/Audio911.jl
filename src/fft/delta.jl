# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
"""
    DeltaSetup <: AbstractSetup

Configuration structure for delta (and delta-delta) coefficient computation.

# Fields
- `sr::Int64`: Sample rate of the source audio signal, in Hz.
- `delta_length::Int64`: Length of the regression window used to compute delta
  coefficients. Must be odd; default is `9`.
- `source::Symbol`: Filtering convention applied to the input matrix:
  - `:standard` — filter along rows (MATLAB-compatible ordering).
  - `:transposed` — filter along columns (AudioFlux-compatible ordering).
"""
struct DeltaSetup <: AbstractSetup
    sr::Int64
    delta_length::Int64
    source::Symbol
end

# ---------------------------------------------------------------------------- #
#                                     mfcc                                     #
# ---------------------------------------------------------------------------- #
"""
    Delta{T} <: AbstractDelta

First-order delta coefficients of a cepstral feature sequence (e.g., MFCCs).

Delta coefficients estimate the **local time derivative** of a feature sequence and
capture the dynamics of the signal over time. They are commonly appended to static
cepstral features (e.g., MFCCs) to enrich the feature vector with temporal
information, improving performance in speech and audio recognition tasks.

Computation uses a **linear regression** over a sliding window of length
`delta_length`, which is equivalent to convolving the feature sequence with a
set of regression weights. The filter is applied causally using `DSP.filt`.

# Fields
- `spec::AbstractArray{T}`: Delta coefficient array, same shape as the input
  feature matrix.
- `info::DeltaSetup`: Configuration metadata (sample rate, window length,
  source convention).

# Constructors

    Delta{S}(x::AbstractArray{T}; sr, delta_length=9, source=:standard)

Low-level constructor. Computes delta coefficients from a raw feature array.

## Arguments
- `x::AbstractArray{T}`: Input feature matrix (e.g., MFCC coefficients).
- `sr::Int64`: Sample rate of the source signal, in Hz.
- `delta_length::Int64=9`: Regression window length. Must be an odd integer ≥ 3.
  Longer windows produce smoother derivatives but introduce more temporal lag.
- `source::Symbol=:standard`: Filtering axis convention:
  - `:standard` — filter applied along rows (MATLAB style).
  - `:transposed` — filter applied along columns (AudioFlux style).

---

    Delta(x::AbstractCepstrum; kwargs...)

Convenience constructor. Computes delta coefficients directly from a cepstral
feature object (e.g., an MFCC struct). Sample rate is extracted automatically.

---

    Delta(x::Delta{S}; kwargs...)

Convenience constructor. Computes delta coefficients from an existing `Delta`
object, yielding **delta-delta** (acceleration) coefficients when used standalone.
For the standard paired workflow, prefer [`DeltaDelta`](@ref).

# Algorithm
The regression weights ``b`` for a window of half-length ``m = \\lfloor L/2 \\rfloor``
are:

```math
b[k] = \\frac{m - k + 1}{\\sum_{i=1}^{m} i^2}, \\quad k = 1, \\ldots, L
```

The delta of feature sequence ``c_t`` at frame ``t`` is approximated as:

```math
\\Delta c_t \\approx \\frac{\\sum_{k=-m}^{m} k \\cdot c_{t+k}}{\\sum_{k=1}^{m} k^2}
```

This is equivalent to fitting a straight line through the `delta_length` neighboring
frames and taking its slope, providing a noise-robust finite-difference estimate.
"""
struct Delta{T} <: AbstractDelta
    spec::AbstractArray{T}
    info::DeltaSetup

    function Delta{S}(
        x::AbstractArray{T};
        sr::Int64,
        delta_length::Int64=9,
        source::Symbol=:standard,
    ) where {T<:AbstractFloat, S}
        # validate delta_length
        delta_length > 2 || throw(ArgumentError("delta_length must be > 2, got $delta_length."))
        isodd(delta_length)  || throw(ArgumentError("delta_length must be odd, got $delta_length."))
        
        # define window shape
        m = delta_length ÷ 2
        b = collect(m:-1:(-m)) ./ sum((1:m) .^ 2)

        spec = if source == :transposed
            DSP.filt(b, 1.0, x')'   #:audioflux setting
        else
            DSP.filt(b, 1.0, x)     #:matlab setting
        end
        info = DeltaSetup(sr, delta_length, source)

        new{T}(spec, info)
    end

    Delta(x::AbstractCepstrum; kwargs...) =
        Delta{typeof(x)}(get_data(x); sr=get_sr(x), kwargs...)

    Delta(x::Delta{S}; kwargs...) where {S} =
        Delta{S}(get_data(x); sr=get_sr(x), kwargs...)
end

"""
    DeltaDelta(x::AbstractCepstrum, kwargs...) -> (Delta, Delta)

Compute first-order (delta) and second-order (delta-delta) coefficients from a
cepstral feature sequence in a single call.

Delta-delta coefficients capture the **acceleration** of feature trajectories over
time. Together with static features and first-order deltas, they form the standard
39-dimensional MFCC feature vector widely used in speech and audio recognition
(13 static + 13 Δ + 13 ΔΔ).

# Arguments
- `x::AbstractCepstrum`: Input cepstral feature object (e.g., MFCCs).
- `kwargs...`: Additional keyword arguments forwarded to the [`Delta`](@ref)
  constructor (e.g., `delta_length`, `source`).

# Returns
A tuple `(d1, d2)` where:
- `d1::Delta`: First-order delta coefficients (velocity).
- `d2::Delta`: Second-order delta-delta coefficients (acceleration), computed
  by applying the same delta filter to `d1`.

# See also
[`Delta`](@ref)
"""
function DeltaDelta(x::AbstractCepstrum, kwargs...)
    d1 = Delta{typeof(x)}(get_data(x); sr=get_sr(x), kwargs...)
    d2 = Delta{typeof(d1)}(get_data(d1); sr=get_sr(x), kwargs...)
    return d1, d2
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
"""
    Base.eltype(::Delta{T}) -> Type

Return the element type of the delta spectrogram data.
"""
Base.eltype(::Delta{T}) where {T} = T

"""
    get_data(m::Delta) -> Matrix

Get the delta spectrogram data matrix, transposed to (nframes × nbands).
"""
@inline get_data(m::Delta)  = m.spec

"""
    get_setup(m::Delta) -> MelSpecSetup

Get the configuration metadata for the delta spectrogram.
"""
@inline get_setup(m::Delta) = m.info

"""
    get_sr(m::Delta) -> Int64

Return the sample rate (in Hz) associated with the Delta spectrogram.
"""
@inline get_sr(m::Delta) = m.info.sr
