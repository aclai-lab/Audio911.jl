# include("windows.jl")

function buffer(
        x::AbstractVector{Float64},
        win_length::Int64,
        hop_length::Int64
)
    x_length = size(x, 1)
    num_hops = floor(Int, (x_length - win_length) / hop_length) + 1

    y = zeros(Float64, win_length, num_hops)

    for j in 1:num_hops
        for i in 1:win_length
            y[i, j] = x[i + hop_length * (j - 1)]
        end
    end

    return y
end # function buffer

function logEnergyCoeffs(
        x::AbstractArray{T}
) where {T <: Real}
    DT = eltype(x)
    E = sum(x .^ 2, dims = 1) # eleva tutti gli elementi ^2 e li somma per colonna
    E[E .== 0] .= floatmin(DT) # se un valore è zero, lo sostituisce col valore più piccolo positivo possibile, in accordo col tipo utilizzato
    logE = log.(E) # fa il log di tutti gli elementi
end # logEnergyCoeffs

function windowing(
        x::Union{AbstractVector{T}, AbstractArray{T}},
        fftLength::Int64 = 256,
        winType::Symbol = :hann,
        winParam::Symbol = :symmetric,
        logEnergy::Bool = false
) where {T <: Real}
    xLength = size(x, 1) # lunghezza audio
    nChan = size(x, 2) # numero canali (mono, stereo)
    DT = eltype(x) # restituisce il tipo degli elementi

    # parto con un if then ma sarebbe bello implementare un Dict
    if (winType == :hann || winType == :hamming || winType == :blackman ||
        winType == :flattopwin)
        win, hL = gencoswin(winType, fftLength, winParam)
    elseif (winType == :rect)
        win, hL = rectwin(fftLength)
    end

    wincast = convert.(DT, win[:]) # casta al tipo originale del file audio
    nHops = Integer(floor((xLength - fftLength) / hL) + 1) # numero delle window necessarie
    y = buffer(x, fftLength, hL)

    # calcola il vettore coefficienti energia. in matlab lo fa qui, ma forse è logico spostarlo dopo il windowing?
    if (logEnergy == true)
        logE = logEnergyCoeffs(y)
    end

    return y' .* wincast'
end # function windowing

function fade(
        x::Union{AbstractVector{T}, AbstractArray{T}},
        fftLength::Int64,
        type::Symbol
) where {T <: Real}
    xLength = size(x, 1) # lunghezza audio
    nChan = size(x, 2) # numero canali (mono, stereo)
    DT = eltype(x) # restituisce il tipo degli elementi

    win, hL = gencoswin(:hann, fftLength, :symmetric)

    wincast = convert.(DT, win[:]) # casta al tipo originale del file audio

    # for c = 1:numChan
    if (type == :in)
        for w in 1:Int(round(fftLength / 2))
            x[w] *= wincast[w]
        end
    elseif (type == :out)
        for w in Int(round(fftLength / 2)):-1:fftLength
            x[end - fftLength + w] *= wincast[w]
        end
    end
    # end

    return x
end # function fade

#------------------------------------------------------------------------------#
#                                   windowing                                  #
#------------------------------------------------------------------------------#
function _get_frames(
	x::AbstractVector{Float64};
	win_type::Tuple{Symbol, Symbol},
	win_length::Int64,
	overlap_length::Int64,
)
    hop_length = win_length - overlap_length
	frames = buffer(x, win_length, hop_length)
	win, _ = gencoswin(win_type[1], win_length, win_type[2])

    n_col = div(size(x, 1) - overlap_length, hop_length)
    col_offsets = (0:(n_col-1)) .* hop_length

	return frames, win, frames .* win, n_col, col_offsets
end

function _get_frames(x::AbstractVector{Float64}, s::StftSetup)
	_get_frames(x, win_type = s.win_type, win_length = s.win_length, overlap_length = s.overlap_length)
end

function _get_frames(x::AbstractVector{Float64}, a::AudioSetupDev)
	_get_frames(x, a.stft)
end

function get_frames(
	x::AbstractVector{<:AbstractFloat},
	win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
	win_length::Int64 = 256,
	overlap_length::Int64 = 128,
)
	@assert 0 < overlap_length < win_length "Overlap length must be < window length."

	frames, win, wframes, n_col, col_offsets = _get_frames(
		eltype(x) == Float64 ? x : Float64.(x),
		win_type = win_type,
		win_length = win_length,
		overlap_length = overlap_length,
	)
end

get_frames!(a::AudioObj; kwargs...) = a.data.frames, a.setup.window = get_frames(a.data.x; kwargs...)