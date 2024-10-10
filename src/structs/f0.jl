# ---------------------------------------------------------------------------- #
#                             fundamental frequency                            #
# ---------------------------------------------------------------------------- #
struct F0Setup
	method::Symbol
	freqrange::Tuple{Int, Int}
	mflength::Int
end

struct F0Data
    f0::AbstractVector{<:AbstractFloat}
end

struct F0
    sr::Int
    setup::F0Setup
    data::F0Data
end

# ---------------------------------------------------------------------------- #
#                                    utilities                                 #
# ---------------------------------------------------------------------------- #
function get_candidates(domain::AbstractArray{T}, freqrange::NTuple{2,Int}) where T<:AbstractFloat
    range = freqrange[1]:freqrange[2]
    peaks, locs = findmax(view(domain, range, :), dims=1)
	locs = freqrange[1] .+ map(i -> i[1], locs) .- 1
    return vec(peaks), locs
end


function i_clip(x::AbstractVector{<:AbstractFloat}, range::Tuple{Int, Int})
	x[x.<range[1]] .= range[1]
	x[x.>range[2]] .= range[2]

	return x
end

function compute_residual(x::AbstractArray{T}, invrt::Vector{Vector{T}}, nwin::Int, nhop::Int, hops::Int) where T <: AbstractFloat
    x_length = length(x)
    
	residual = zeros(T, x_length)
    @inbounds for i in 1:hops
        idx1 = 1 + nhop * (i - 1)
        idx2 = min(nwin + nhop * (i - 1), x_length)
        @views residual[idx1:idx2] .+= invrt[i][1:min(idx2-idx1+1, nwin)]
    end

	residual_buff = (let
		idx1 = 1 + nhop * (hop - 1)
		idx2 = min(nwin + nhop * (hop - 1), x_length)
		temp = @view residual[idx1:idx2]
		[temp; zeros(T, max(0, nwin - length(temp)))]
	end for hop in 1:hops)

	win = _get_window(:blackman, nwin, true)
	_get_wframes(residual_buff, win)
end

# ---------------------------------------------------------------------------- #
#                             fundamental frequency                            #
# ---------------------------------------------------------------------------- #
function _get_f0(
	x::AbstractArray{<:AbstractFloat},
	sr::Int;
	x_length::Int,
	nwin::Int,
	noverlap::Int,
	method::Symbol = :nfc,
	freqrange::Tuple{Int, Int} = (50, 400),
	mflength::Int = 1,
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

	elseif method == :srh
		nhop = nwin - noverlap
		hops = ceil(Int, (x_length - nwin) / nhop)

		order = 12
		a = combinedims(collect(vcat(1.0, lpc(col, order-1, LPCLevinson())[1]) for col in eachcol(x)))'
		invrt = collect(filt(row, [1.0], x[:, i]) for (i, row) in enumerate(eachrow(a)))

		residual_buff = compute_residual(x, invrt, nwin, nhop, hops)

		maxbin = 5 * freqrange[2]
		residual_pad = (vcat(r, zeros(eltype(r), sr - nwin)) for r in residual_buff)
		res = fft(combinedims(collect(residual_pad)), 1)

		e = @views abs.(res[1:maxbin, :])

		domain = zeros(eltype(e), freqrange[2], size(e, 2))
		freqr = freqrange[1]:freqrange[2]

		@views @. domain[freqr, :] = 
			e[freqr, :] + 
			e[2*freqr, :] - e[round(Int, 1.5*freqr), :] + 
			e[3*freqr, :] - e[round(Int, 2.5*freqr), :] + 
			e[4*freqr, :] - e[round(Int, 3.5*freqr), :] + 
			e[5*freqr, :] - e[round(Int, 4.5*freqr), :]

		_, locs = get_candidates(domain, freqrange)
		f0 = vec(Float64.(locs))

	# elseif method == :pef
	# elseif method == :cep
	# elseif method == :lhs

	end

	F0(sr, F0Setup(method, freqrange, mflength), F0Data(f0))
end

function Base.show(io::IO, f0::F0)
    println(io, "F0 Estimation:")
    println(io, "  Sample Rate: $(f0.sr) Hz")
    println(io, "  Method: $(f0.setup.method)")
    println(io, "  Frequency Range: $(f0.setup.freqrange) Hz")
    println(io, "  F0: $(length(f0.data.f0)) points")
end

function Plots.plot(f0::F0)
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

get_f0(source::Stft; kwargs...) = _get_f0(combinedims(collect(source.data.frames)), source.sr; x_length=source.x_length, nwin=source.setup.nwin, noverlap=source.setup.noverlap, kwargs...)
