# using Statistics

# include("../windowing/windows.jl")
# include("../windowing/windowing.jl")

function signal_to_rms(
	x::Union{AbstractVector{T}, AbstractArray{T}},
) where {T <: AbstractFloat}
	# calcola power e rms
	sqrt.(mean(abs.(x) .^ 2, dims = 2))
end

function signal_to_db(
	x::Union{AbstractVector{T}, AbstractArray{T}},
) where {T <: AbstractFloat}
	magnitude = abs.(x)
	power = magnitude .^ 2

	ref_value = maximum(magnitude)^2
	amin = 1e-5^2

	# power to db
	db_log = 10.0 * log10.(max.(maximum.(power), amin))
	db_log .-= 10.0 * log10(max(amin, ref_value))
	return db_log
end

function trim_audio(
	x::AbstractVector{Float64};
	fft_length::Int64 = 1024,
	threshold::Real = 60,
)

	# crea la matrice delle finestre: le righe sono le finestre.
	y = windowing(x, fft_length, :rect, :symmetric, false)

	# converte l'audio in rms e poi in db
	db = signal_to_db(signal_to_rms(y))
	framesValidated = db .> -threshold

	signal = zeros(eltype(x), 0)
	silence = zeros(eltype(x), 0)
	flag = :start

	for i in eachindex(framesValidated)
		# case 1: starting signal portion
		if (framesValidated[i] && (flag == :silence || flag == :start))
			fade_y = fade(y[i, :], fft_length, :in)
			signal = vcat(signal, y[i, :])

			if (flag == :silence)
				silence = fade(silence, fft_length, :out)
			end

			flag = :signal

			# case 2: continuing signal portion
		elseif (framesValidated[i] && flag == :signal)
			signal = vcat(signal, y[i, :])

			# case 3: starting silence portion
		elseif (!framesValidated[i] && (flag == :signal || flag == :start))
			fade_y = fade(y[i, :], fft_length, :in)
			silence = vcat(silence, y[i, :])

			if (flag == :signal)
				signal = fade(signal, fft_length, :out)
			end

			flag = :silence

			# case 4: continuing silence portion
		else
			silence = vcat(silence, y[i, :])
		end

		# ending with fade out
		if (flag == :silence)
			silence = fade(silence, fft_length, :out)
		else
			signal = fade(signal, fft_length, :out)
		end
	end

	signal, silence
end

function trim_audio(
	x::AbstractVector{T};
	fft_length::Int64 = 1024,
	threshold::Real = 60,
) where {T <: AbstractFloat}
	x_type = eltype(x)
	signal, silence = trim_audio(Float64.(x), fft_length = fft_length, threshold = threshold)
	# x_type.(signal), x_type.(silence)
end

function normalize_audio(x::AbstractVector{T}) where {T <: AbstractFloat}
	if length(x) != 0
		x ./ maximum(abs.(x))
	end
end
