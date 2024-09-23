# ---------------------------------------------------------------------------- #
#                             fundamental frequency                            #
# ---------------------------------------------------------------------------- #
struct F0Setup
	method::Symbol
	freqrange::Tuple{Int64, Int64}
	mflength::Int64
end

struct F0Data
    f0::AbstractVector{<:AbstractFloat}
end

struct F0
    sr::Int64
    setup::F0Setup
    data::F0Data
end

# ---------------------------------------------------------------------------- #
#                                    utilities                                 #
# ---------------------------------------------------------------------------- #
function get_candidates(domain::AbstractArray{Float64}, edge::Tuple{Int64, Int64})
	range = collect(edge[1]:edge[end])
	lower = edge[1]

	peaks, locs = findmax(domain[range, :], dims = 1)
	locs = lower .+ map(i -> i[1], locs) .- 1

	return peaks, locs
end

function i_clip(x::AbstractVector{Float64}, range::Tuple{Int64, Int64})
	x[x.<range[1]] .= range[1]
	x[x.>range[2]] .= range[2]

	return x
end

# ---------------------------------------------------------------------------- #
#                             fundamental frequency                            #
# ---------------------------------------------------------------------------- #
function _get_f0(;
	x::AbstractArray{Float64},
	sr::Int64,
	nwin::Int64,
	method::Symbol = :nfc,
	freqrange::Tuple{Int64, Int64} = (50, 400),
	mflength::Int64 = 1,
)
	if method == :nfc
		edge = (round.(Int, sr ./ reverse(freqrange)))
		r = size(x, 1)
		mxl = min(edge[2], r - 1)
		m2 = nextpow(2, 2 * nwin - 1)
		y_m2 = zeros(m2, size(x, 2))
		for i in axes(x, 1)
			y_m2[i, :] = x[i, :]
		end
		c1 = real.(ifft(abs.(fft(y_m2, (1,))) .^ 2, (1,))) ./ sqrt(m2)

		Rt  = [c1[m2-mxl.+collect(1:mxl), :]; c1[1:mxl+1, :]]
		lag = Rt[(edge[2] + 1 + edge[1]):end, :]

		# energy of original signal, y
		yRMS = sqrt.(Rt[edge[2] + 1, :])

		# normalize lag domain energy by input signal energy
		lag = lag ./ yRMS'

		# zero-pad domain so time locs can be easily interpreted.
		domain = [zeros(edge[1] - 1, size(lag, 2)); lag]

		# peak picking
		_, locs = get_candidates(domain, edge)

		# convert lag domain to frequency
		f0 = vec(sr ./ locs)

		## TODO
		# elseif f0.method == :srh
		# elseif f0.method == :pef
		# elseif f0.method == :cep
		# elseif f0.method == :lhs

		F0(sr, F0Setup(method, freqrange, mflength), F0Data(f0))
	end
end

function Base.show(io::IO, f0::F0)
    println(io, "F0 Estimation:")
    println(io, "  Sample Rate: $(f0.sr) Hz")
    println(io, "  Method: $(f0.setup.method)")
    println(io, "  Frequency Range: $(f0.setup.freqrange) Hz")
    println(io, "  F0: $(length(f0.data.f0)) points")
end

function Base.display(f0::F0)
    time = (0:length(f0.data.f0)-1)
    plot(time, f0.data.f0, 
         title="Fundamental Frequency Estimation",
         xlabel="Frame",
         ylabel="Frequency (Hz)",
         label="F0",
         linewidth=2,
         ylim=(0, maximum(f0.data.f0) * 1.1),
         legend=:none)
end

function get_f0(;
	source::Stft,
	kwargs...
)
	_get_f0(; x=source.data.frames, sr=source.sr, nwin=source.setup.nwin, kwargs...)
end