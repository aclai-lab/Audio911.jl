# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
"""
    AbstractWinFunction

Base type for windowing function implementations.
"""
abstract type AbstractWinFunction end

# ---------------------------------------------------------------------------- #
#                                win functions                                 #
# ---------------------------------------------------------------------------- #
"""
    WinFunction <: AbstractWinFunction

Callable wrapper for windowing algorithms with parameters.

# Fields
- `func::Function`: The windowing implementation function
- `params::NamedTuple`: Algorithm-specific parameters
"""
struct WinFunction <: AbstractWinFunction
    func   :: Function
    params :: NamedTuple
end
# Make it callable - npoints is passed at execution time
(w::WinFunction)(npoints::Int; kwargs...) = w.func(npoints; w.params..., kwargs...)

"""
    MovingWindow(; window_size::Int, window_step::Int) -> WinFunction

Create a moving window that slides across the time series.

# Parameters
- `window_size`: Number of time points in each window
- `window_step`: Step size between consecutive windows

# Example
    win = MovingWindow(window_size=10, window_step=5)
    intervals = win(100)  # For 100-point time series
"""
function MovingWindow(;
    window_size::Int,
    window_step::Int,
)
    WinFunction(movingwindow, (;window_size, window_step))
end

"""
    WholeWindow() -> WinFunction

Create a single window encompassing the entire time series.
Useful for global feature extraction without temporal partitioning.

# Example
    win = WholeWindow()
    intervals = win(100)  # Returns [1:100]
"""
WholeWindow(;) = WinFunction(wholewindow, (;))

"""
    SplitWindow(; nwindows::Int) -> WinFunction

Divide the time series into equal non-overlapping segments.

# Parameters
- `nwindows`: Number of equal-sized windows to create

# Example
    win = SplitWindow(nwindows=4)
    intervals = win(100)  # Four 25-point windows
"""
SplitWindow(;nwindows::Int) = WinFunction(splitwindow, (; nwindows))

"""
    AdaptiveWindow(; nwindows::Int, relative_overlap::AbstractFloat) -> WinFunction

Create overlapping windows with adaptive sizing based on series length.

# Parameters
- `nwindows`: Target number of windows
- `relative_overlap`: Fraction of overlap between adjacent windows (0.0-1.0)

# Example
    win = AdaptiveWindow(nwindows=3, relative_overlap=0.1)
    intervals = win(100)  # Three adaptive windows with 10% overlap
"""
function AdaptiveWindow(;
    nwindows::Int,
    relative_overlap::AbstractFloat,
)
    WinFunction(adaptivewindow, (; nwindows, relative_overlap))
end



# function logEnergyCoeffs(
#         x::AbstractArray{T}
# ) where {T <: Real}
#     DT = eltype(x)
#     E = sum(x .^ 2, dims = 1) # eleva tutti gli elementi ^2 e li somma per colonna
#     E[E .== 0] .= floatmin(DT) # se un valore è zero, lo sostituisce col valore più piccolo positivo possibile, in accordo col tipo utilizzato
#     logE = log.(E) # fa il log di tutti gli elementi
# end # logEnergyCoeffs

# function windowing(
#         x::Union{AbstractVector{T}, AbstractArray{T}},
#         fftLength::Int64 = 256,
#         winType::Symbol = :hann,
#         winParam::Symbol = :symmetric,
#         logEnergy::Bool = false
# ) where {T <: Real}
#     xLength = size(x, 1) # lunghezza audio
#     nChan = size(x, 2) # numero canali (mono, stereo)
#     DT = eltype(x) # restituisce il tipo degli elementi

#     # parto con un if then ma sarebbe bello implementare un Dict
#     if (winType == :hann || winType == :hamming || winType == :blackman ||
#         winType == :flattopwin)
#         win, hL = gencoswin(winType, fftLength, winParam)
#     elseif (winType == :rect)
#         win, hL = rectwin(fftLength)
#     end

#     wincast = convert.(DT, win[:]) # casta al tipo originale del file audio
#     nHops = Integer(floor((xLength - fftLength) / hL) + 1) # numero delle window necessarie
#     y = buffer(x, fftLength, hL)

#     # calcola il vettore coefficienti energia. in matlab lo fa qui, ma forse è logico spostarlo dopo il windowing?
#     if (logEnergy == true)
#         logE = logEnergyCoeffs(y)
#     end

#     return y' .* wincast'
# end # function windowing

# function fade(
#         x::Union{AbstractVector{T}, AbstractArray{T}},
#         fftLength::Int64,
#         type::Symbol
# ) where {T <: Real}
#     xLength = size(x, 1) # lunghezza audio
#     nChan = size(x, 2) # numero canali (mono, stereo)
#     DT = eltype(x) # restituisce il tipo degli elementi

#     win, hL = gencoswin(:hann, fftLength, :symmetric)

#     wincast = convert.(DT, win[:]) # casta al tipo originale del file audio

#     # for c = 1:numChan
#     if (type == :in)
#         for w in 1:Int(round(fftLength / 2))
#             x[w] *= wincast[w]
#         end
#     elseif (type == :out)
#         for w in Int(round(fftLength / 2)):-1:fftLength
#             x[end - fftLength + w] *= wincast[w]
#         end
#     end
#     # end

#     return x
# end # function fade

#------------------------------------------------------------------------------#
#                                   windowing                                  #
#------------------------------------------------------------------------------#
# function _get_frames(
# 	x::AbstractVector{Float64};
# 	window_type::Tuple{Symbol, Symbol},
# 	window_length::Int64,
# 	overlap_length::Int64,
# )
# 	frames = buffer(x, window_length, window_length - overlap_length)
# 	window, _ = gencoswin(window_type[1], window_length, window_type[2])

# 	return frames, window
# end

# function _get_frames(x::AbstractVector{Float64}, s::AudioSetup)
# 	_get_frames(x, window_type = s.window_type, window_length = s.window_length, overlap_length = s.overlap_length)
# end

# function get_frames(
# 	x::AbstractVector{<:AbstractFloat},
# 	window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# 	window_length::Int64 = 256,
# 	overlap_length::Int64 = 128,
# )
# 	@assert 0 < overlap_length < window_length "Overlap length must be < window length."

# 	frames, window = _get_frames(
# 		eltype(x) == Float64 ? x : Float64.(x),
# 		window_type = window_type,
# 		window_length = window_length,
# 		overlap_length = overlap_length,
# 	)
# end

# get_frames!(a::AudioObj; kwargs...) = a.data.frames, a.setup.window = get_frames(a.data.x; kwargs...)

#------------------------------------------------------------------------------#
#                                    frames                                    #
#------------------------------------------------------------------------------#
function get_frames(
	X    :: AudioFile;
    win  :: WinFunction=AdaptiveWindow(nwindows=3, relative_overlap=0.1),
	type :: Tuple{Symbol, Symbol}=(:hann, :periodic),
)
    # run the windowing algo and set windows indexes
    intervals = win(length(X))

    map(intervals) do interval
        frame      = data(X)[interval]
        window, _  = gencoswin(type[1], length(frame), type[2])
        frame .* window
    end
end

# function _get_frames2(x::AbstractVector{Float64}, s::AudioSetup)
# 	_get_frames2(x, window_type = s.window_type, window_length = s.window_length, overlap_length = s.overlap_length)
# end

# function get_frames2(
# 	x::AbstractVector{<:AbstractFloat},
# 	window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# 	window_length::Int64 = 256,
# 	overlap_length::Int64 = 128,
# )
# 	@assert 0 < overlap_length < window_length "Overlap length must be < window length."

# 	frames = _get_frames2(
# 		eltype(x) == Float64 ? x : Float64.(x),
# 		window_type = window_type,
# 		window_length = window_length,
# 		overlap_length = overlap_length,
# 	)
# end

# get_frames2!(a::AudioObj; kwargs...) = a.data.frames, a.setup.window = get_frames2(a.data.x; kwargs...)
