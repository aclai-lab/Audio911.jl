# using FFTW

function get_candidates(
	domain::AbstractArray{T},
	edge::AbstractVector{Int64},
) where {T <: AbstractFloat}

	range = collect(edge[1]:edge[end])
	lower = edge[1]

	peaks, locs = findmax(domain[range, :], dims = 1)
	locs = lower .+ map(i -> i[1], locs) .- 1

	return peaks, locs
end

function i_clip(
	x::AbstractVector{T},
	range::Vector{Int64} = [50, 400],
) where {T <: AbstractFloat}

	x[x.<range[1]] .= range[1]
	x[x.>range[2]] .= range[2]

	return x
end

function f0(
	data::signal_data,
	setup::signal_setup,
	method::Symbol = :nfc,
	f0_range::Vector{Int64} = [50, 400],
	median_filter_length::Int64 = 1,
)
	len_x = size(data.x, 1)
	hoplength = setup.window_length - setup.overlap_length
	num_hops_final = Int(floor((len_x - setup.window_length) / hoplength)) + 1

	if method == :srh
		N = Int(round(0.025 * setup.sr))
		hop_size = Int(round(0.005 * setup.sr))
	else
		N = setup.window_length
		hop_size = hoplength
	end
	num_hops = Int(ceil((len_x - N) / hop_size)) + 1

	num_to_pad = (N + num_hops * hop_size) - len_x - 1

	# split in windows
	y = buffer([data.x; zeros(num_to_pad)], N, hop_size)
	extra_params = Dict(:num_candidates => 1, :min_peak_distance => 1)

	if method == :nfc
		edge = Int.(round.(setup.sr ./ reverse(f0_range)))
		r = size(y, 1)
		mxl = min(edge[end], r - 1)
		m2 = nextpow(2, 2 * r - 1)
		y_m2 = zeros(m2, size(y, 2))
		for i in axes(y, 1)
			y_m2[i, :] = y[i, :]
		end
		c1 = real.(ifft(abs.(fft(y_m2, (1,))) .^ 2, (1,))) ./ sqrt(m2)

		Rt  = [c1[m2-mxl.+collect(1:mxl), :]; c1[1:mxl+1, :]]
		lag = Rt[(edge[end]+1+edge[1]):end, :]

		# energy of original signal, y
		yRMS = sqrt.(Rt[edge[end]+1, :])

		# normalize lag domain energy by input signal energy
		lag = lag ./ yRMS'

		# zero-pad domain so time locs can be easily interpreted.
		domain = [zeros(edge[1] - 1, size(lag, 2)); lag]

		# peak picking
		peaks, locs = get_candidates(domain, edge)

		# convert lag domain to frequency
		frq_0 = setup.sr ./ locs
		# conf = peak./sum(abs(peak),2)

		## TODO
		# elseif method == :srh
		# elseif method == :pef
		# elseif method == :cep
		# elseif method == :lhs
	end

	# force pitch estimate inside band edges
	frq_0 = i_clip(vec(frq_0), f0_range)
	data.f0 = vec(frq_0[1:num_hops_final, :])
end