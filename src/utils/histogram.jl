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
	xmin, xmax = extrema(filter(x -> !isnan(x), feature))
	range_feature = xmax - xmin
	raw_bins_width = range_feature / hist_bins

	edges = binpicker(xmin, xmax, hist_bins, raw_bins_width)

	n, _ = histcountindices(feature, edges)

	return n, edges
end