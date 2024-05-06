function moving_mean(
	x::Vector{Float64},
	w::Int64,
)
	# w must be odd!
	x_length = size(x, 1)
	m = zeros(x_length)
	pad = Int(floor(w / 2))
	for i in eachindex(x)
		index1 = i - pad < 1 ? 1 : i - pad
		index2 = i + pad > x_length ? x_length : i + pad
		m[i] = median(x[index1:index2])
	end
	return m
end

function binpicker(
	xmin::Float64,
	xmax::Float64,
	nbins::Int64,
	raw_bins_width::Float64,
)
	xscale = max(abs(xmin), abs(xmax))
	xrange = xmax - xmin
	raw_bins_width = max(raw_bins_width, eps(xscale))

	# if the data are not constant, place the bins at "nice" locations
	if xrange > max(sqrt(eps(Float64)) * xscale, floatmin(Float64))
		# choose the bin width as a "nice" value.
		pow_of_ten = 10 .^ floor(log10(raw_bins_width))   # next lower power of 10
		rel_size = raw_bins_width / pow_of_ten           # guaranteed in [1, 10)

		# automatic rule specified
		# if isempty(nbins)      nbins viene passato quindi questo codice non serve
		#     if  rel_size < 1.5
		#         bin_width = 1*pow_of_ten
		#     elseif rel_size < 2.5
		#         bin_width = 2*pow_of_ten
		#     elseif rel_size < 4
		#         bin_width = 3*pow_of_ten
		#     elseif rel_size < 7.5
		#         bin_width = 5*pow_of_ten
		#     else
		#         bin_width = 10*pow_of_ten
		#     end

		# put the bin edges at multiples of the bin width, covering x.
		# the actual number of bins used may not be exactly equal to the requested rule.
		# left_edge = max(min(bin_width*floor(xmin ./ bin_width), xmin),-realmax(class(xmax)))
		# nbinsActual = max(1, ceil((xmax-left_edge) ./ bin_width))
		# rightEdge = min(max(left_edge + nbinsActual.*bin_width, xmax),realmax(class(xmax)))

		# number of bins specified
		# else
		# temporarily set a raw bin_width to a nice power of 10.
		# bin_width will be set again to a different value if more than 1 bin.
		bin_width = pow_of_ten * floor(rel_size)
		# set the left edge at multiples of the raw bin width.
		# then adjust bin width such that all bins are of the same size and xmax fall into the rightmost bin.
		left_edge = max(min(bin_width * floor(xmin ./ bin_width), xmin), -floatmax(Float64))
		if nbins > 1
			ll = (xmax - left_edge) / nbins      # bin_width lower bound, xmax
			# on right edge of last bin
			ul = (xmax - left_edge) / (nbins - 1)  # bin_width upper bound,
			# xmax on left edge of last bin
			p10 = 10^floor(log10(ul - ll))
			bin_width = p10 * ceil(ll ./ p10)    # bin_width-ll < p10 <= ul-ll
			# thus, bin_width < ul
		end

		nbins_actual = nbins
		right_edge = min(
			max(left_edge + nbins_actual .* bin_width, xmax), floatmax(Float64))
		# end

	else # the data are nearly constant
		# for automatic rules, use a single bin.
		if isempty(nbins)
			nbins = 1
		end

		# there's no way to know what scale the caller has in mind, just create
		# something simple that covers the data.

		# make the bins cover a unit width, or as small an integer width as
		# possible without the individual bin width being zero relative to xscale.
		# put the left edge on an integer or half integer below xmin, with the data in the middle 50% of the bin.
		# put the right edge similarly above xmax.
		bin_range = max(1, ceil(nbins * eps(xscale)))
		left_edge = floor(2 * (xmin - bin_range ./ 4)) / 2
		right_edge = ceil(2 * (xmax + bin_range ./ 4)) / 2

		bin_width = (right_edge - left_edge) ./ nbins
		nbins_actual = nbins
	end

	if !isfinite(bin_width)
		# if binWidth overflows, don't worry about nice bin edges anymore
		edges = LinRange(left_edge, right_edge, nbins_actual + 1)
	else
		edges = union(
			left_edge, left_edge .+ (1:(nbins_actual-1)) .* bin_width, right_edge)
		step = round(minimum(diff(edges)), digits = 8)
		edges = range(edges[1], edges[end], step = step)
	end

	return edges
end

function histcounts(
	feature::Vector{Float64},
	hist_bins::Int64,
)
	edgestransposed = false

	xmin = minimum(feature)
	xmax = maximum(feature)
	range_feature = xmax - xmin
	raw_bins_width = range_feature / hist_bins

	edges = binpicker(xmin, xmax, hist_bins, raw_bins_width)

	n, _ = histcountindices(feature, edges)

	return n, edges
end

function f_peaks(n::Vector{T}) where {T <: AbstractFloat}
	z7 = zeros(Float64, 7)
	z8 = zeros(Float64, 8)
	z3 = zeros(Float64, 3)

	nn = vcat(n, z7, n, z7, n, z8, n, z7, n, z7, n)

	n[end] = 0
	temp = repeat([z3; n; z3], 6, 1)

	b = all(reshape(nn .< temp, (Int(round(length(nn) / 6)), 6)), dims = 2)

	peaks_idx = []
	for i in eachindex(b)
		if b[i]
			push!(peaks_idx, i)
		end
	end
	peaks_idx = peaks_idx .- 3

	return peaks_idx
end

function get_threshs_from_feature(
	feature::Vector{Float64},
	bins::Int64,
	type::Symbol,
)
	# get histogram
	hist_bins = round(Int, length(feature) / bins)
	# at leat 10 histogram
	hist_bins = max(10, hist_bins)

	m_feature = mean(feature)
	n_feature, edges_feature = histcounts(feature, hist_bins)

	# working with spectral spread
	if type == :specspread

		# delete the first indices of the histogram if the first bin contains 0s. We put them there when applying the threshold.
		if edges_feature[1] == 0
			n_feature = n_feature[2:end]
			edges_feature = edges_feature[2:end]
		end
		# set minimum values
		minval = m_feature / 2
	else
		minval = minimum(feature)
	end

	# find local maxima for energy
	peaks_idx = f_peaks(Float64.(n_feature))

	# working with energy
	if type == :energy
		# check to make sure that no peaks are beyond halfway through the histogram.
		# high maxima values later in the histogram makes very high threshold values.
		if length(peaks_idx) >= 2 && all(maximum(n_feature[peaks_idx[1:2]]) > m_feature)
			# Remove all maxima that are far along
			if edges_feature[peaks_idx[2]] > m_feature
				peaks_idx = peaks_idx[1]
			elseif edges_feature[peaks_idx[1]] > m_feature
				deleteat!(peaks_idx, 1)
			end
		end
	end

	if length(peaks_idx) == 0
		M1 = m_feature / 2
		M2 = minval
	elseif length(peaks_idx) == 1
		eF0 = vcat(collect(edges_feature), 0)
		M1 = 0.5 * (vcat(0, collect(edges_feature)) .- eF0) + eF0
		M1 = M1[peaks_idx.+1]
		M2 = minval
	else
		eF0 = vcat(collect(edges_feature), 0)
		AA = 0.5 * (vcat(0, collect(edges_feature)) .- eF0) + eF0
		M2 = AA[peaks_idx[1]+1]
		M1 = AA[peaks_idx[2]+1]
	end

	return M1, M2
end

function spectral_spread(
	x::Vector{Float64},
	sr::Int64;
	fft_length::Int64,
	window_length::Int64,
	overlap_length::Int64,
	window_norm::Bool = true,
	spectrum_type::Symbol = :magnitude,
)
	X = audio_features_obj(
		x, sr,
		fft_length = fft_length,
		window_length = window_length,
		overlap_length = overlap_length,
		window_norm = window_norm,
		spectrum_type = spectrum_type,
	)
	s = X.get_lin_spec()

	freq = X.data.lin_frequencies

	sum_x1 = vec(sum(s, dims = 1))
	spectral_centroid = vec(sum(s .* freq, dims = 1) ./ sum_x1')
	spectral_centroid = replace!(spectral_centroid, NaN => 0)
	higher_moment_tmp = freq .- spectral_centroid'

	spectral_spread = vec(sqrt.(sum((higher_moment_tmp .^ 2) .* s, dims = 1) ./ sum_x1'))

	return spectral_spread
end

function speech_detector(
	x_in::AbstractVector{Float64},
	sr::Int64;        #thresholds
)
	window, _ = gencoswin(:hann, Int(round(0.03 * sr)), :periodic)
	frame_length = size(window, 1)
	merge_distance = frame_length * 5

	## helper detect speech
	W = 5                              # weight for finding local maxima
	bins = 15                          # number of bins for histograms
	spectral_spread_threshold = 0.05     # threshold used for not counting spectral spread when under this energy
	lower_spread_threshold_factor = 0.8   # factor to lower the spectral spread threshold.
	smoothing_filter_length = 5          # after getting spectral spread and energy data, these features are smoothed by filters of this length

	#----------------------------------------------------------------------------------#
	#      step 1: extract short-term spectral spread and energy from whole signal     #
	#----------------------------------------------------------------------------------#
	sig_max = maximum(abs.(x_in))

	x = deepcopy(x_in)
	# normalize
	if sig_max > 0
		x = x ./ sig_max
	end

	# buffer signal
	frames = buffer(x, frame_length, frame_length)

	# determine short term energy
	energy = vec(window' .^ 2 * frames .^ 2)
	# filter the short term energy twice
	filtered_energy = moving_mean(
		moving_mean(energy, smoothing_filter_length), smoothing_filter_length)
	# get spectral spread
	spec_spread = spectral_spread(
		x,
		sr,
		fft_length = 2 * frame_length,
		window_length = frame_length,
		overlap_length = 0,
		spectrum_type = :magnitude,
	)
	# normalize the feature
	spec_spread = spec_spread / (sr / 2)
	# set spectral spread value to 0 for frames with low energy
	spec_spread[energy.<spectral_spread_threshold] .= 0
	# filter spectral spread twice
	filtered_spread = moving_mean(
		moving_mean(spec_spread, smoothing_filter_length), smoothing_filter_length)

	#----------------------------------------------------------------------------------#
	#                        step 2: determine thresholds                              #
	#----------------------------------------------------------------------------------#

	# if both energy and spectral spread thresholds are defined
	# TODO
	# calculate energy and spectral spread local bin maxima
	e1_M1, e1_M2 = get_threshs_from_feature(filtered_energy, bins, :energy)
	ss_M1, ss_M2 = get_threshs_from_feature(filtered_spread, bins, :specspread)

	# calculate energy and spectral spread thresholds
	ww = 1 / (W + 1)
	sspread_thresh = ww * (W * ss_M2 + ss_M1[1]) * lower_spread_threshold_factor
	energy_Thresh = ww * (W * e1_M2 + e1_M1[1])

	#----------------------------------------------------------------------------------#
	#                       step 3: apply threshold criterion                          #
	#----------------------------------------------------------------------------------#
	speech_mask = (filtered_spread .> sspread_thresh) .& (filtered_energy .> energy_Thresh)

	#----------------------------------------------------------------------------------#
	#       step 4: Merge frames if they are close as described by MergeDistance       #
	#----------------------------------------------------------------------------------#
	unbuff_out = speech_mask
	frame_length_new = frame_length

	# change frames into data points
	a = repeat(unbuff_out', outer = [frame_length_new, 1])
	unbuff_out_mask = [a[:]; falses(length(x) - length(a), 1)]
	difference = diff([unbuff_out_mask; false], dims = 1)

	# find all changes from speech to silence; return index before change
	idx_m1 = findall(difference .== -1)
	idx_m1 = getindex.(idx_m1, 1)

	# find all changes from silence to speech; return index before change
	if unbuff_out[1] == 0
		idx_p1 = findall(difference .== 1)
		idx_p1 = getindex.(idx_p1, 1)
	else
		idx_p1 = findall(difference .== 1)
		idx_p1 = vcat(1, getindex.(idx_p1, 1))
	end

	# find gaps less than merge distance
	if length(idx_p1) > 1
		testmask = idx_p1[2:end] .- idx_m1[1:(length(idx_p1)-1)] .<= merge_distance
	else
		testmask = falses(0, 1)
	end

	if isempty(idx_p1) || isempty(idx_m1)
		# no speech is found
		outidx = []
	else
		# arrange output
		idx_p2 = idx_p1[2:end, :]
		idx_m2 = idx_m1[1:(length(idx_p1)-1), :]
		amask = .!testmask
		outidx = reshape([idx_p1[1]; idx_p2[amask]; idx_m2[amask]; idx_m1[end]], :, 2)
	end

	y = []
	for i in eachrow(outidx)
		y = [y; x_in[i[1]:i[2]]]
	end

	return Float64.(y), outidx
end

function speech_detector(x_in::AbstractVector{<:AbstractFloat}, sr::Int64)
	speech_detector(Float64.(x_in), sr)
end

# references
# https://github.com/linan2/Voice-activity-detection-VAD-paper-and-code
