# ---------------------------------------------------------------------------- #
#                             fundamental frequency                            #
# ---------------------------------------------------------------------------- #
mutable struct F0
	# setup
    sr::Int64
	method::Symbol
	freq_range::Tuple{Int64, Int64}
	mf_length::Int64
	# data
	f0::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function F0(;
	sr,
	method = :nfc,
	freq_range = (50, 400),
	mf_length = 1,
	f0 = nothing
)
	f0 = F0(
        sr,
		method,
        freq_range,
		mf_length,
		f0
	)
	return f0
end

# ---------------------------------------------------------------------------- #
#                                    utilities                                 #
# ---------------------------------------------------------------------------- #
function get_candidates(
	domain::AbstractArray{Float64},
	edge::Tuple{Int64, Int64},
)
	range = collect(edge[1]:edge[end])
	lower = edge[1]

	peaks, locs = findmax(domain[range, :], dims = 1)
	locs = lower .+ map(i -> i[1], locs) .- 1

	return peaks, locs
end

function i_clip(
	x::AbstractVector{Float64},
	range::Tuple{Int64, Int64}
)
	x[x.<range[1]] .= range[1]
	x[x.>range[2]] .= range[2]

	return x
end

# ---------------------------------------------------------------------------- #
#                             fundamental frequency                            #
# ---------------------------------------------------------------------------- #
function _get_f0(;
	x::AbstractArray{Float64},
	win_length::Int64,
	f0::F0
)
	if f0.method == :nfc
		edge = (round.(Int, f0.sr ./ reverse(f0.freq_range)))
		r = size(x, 1)
		mxl = min(edge[2], r - 1)
		m2 = nextpow(2, 2 * win_length - 1)
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
		f0.f0 = vec(f0.sr ./ locs)

		## TODO
		# elseif f0.method == :srh
		# elseif f0.method == :pef
		# elseif f0.method == :cep
		# elseif f0.method == :lhs

		return f0
	end
end

function Base.show(io::IO, f0::F0)
    println(io, "F0 Estimation:")
    println(io, "  Sample Rate: $(f0.sr) Hz")
    println(io, "  Method: $(f0.method)")
    println(io, "  Frequency Range: $(f0.freq_range) Hz")
    println(io, "  Median Filter Length: $(f0.mf_length)")
    if isnothing(f0.f0)
        println(io, "  F0: Not computed")
    else
        println(io, "  F0: $(length(f0.f0)) points")
    end
end

function Base.display(f0::F0)
    if isnothing(f0.f0)
        println("F0 not computed yet.")
        return
    end

    time = (0:length(f0.f0)-1)

    plot(time, f0.f0, 
         title="Fundamental Frequency Estimation",
         xlabel="Frame",
         ylabel="Frequency (Hz)",
         label="F0",
         linewidth=2,
         ylim=(0, maximum(f0.f0) * 1.1),
         legend=:none)
end

function get_f0(;
	source::Stft,
	kwargs...
)
	_get_f0(x=source.data.frames, win_length=source.setup.win_length, f0=F0(; sr=source.sr, kwargs...))
end