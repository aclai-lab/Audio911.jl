"""
histcounts

input arguments:
x::AbstractVector{Float64};
nbins::Union{Int64, Nothing} = nothing,
norm::Symbol = :none,
				:countdensity, :cumcount, :probability, :percentage, :pdf. :cdf
allow_nan = :false
"""
function binpicker(
	xmin::Float64,
	xmax::Float64,
	nbins::Union{Int64, Nothing},
	binwidth::Union{Int64, Float64},
)
	get_pow10 = (x) -> 10^floor(log10(x))

	xscale = max(abs(xmin), abs(xmax))
	binwidth = max(binwidth, eps(xscale))

	if xmax - xmin > max(sqrt(eps(Float64)) * xscale, Float64(-Inf))
		pow10 = get_pow10(binwidth)
		relsize = binwidth / pow10

		binwidth = if isnothing(nbins)
			conditions = Dict(1.5 => 1, 2.5 => 2, 4 => 3, 7.5 => 5, Float64(Inf) => 10)
			index = findfirst(x -> relsize < x, collect(keys(conditions)))
			conditions[collect(keys(conditions))[index]] * pow10
		else
			pow10 * floor(relsize)
		end

		left_edge = max(min(binwidth * floor(xmin / binwidth), xmin), Float64(-Inf))
		if isnothing(nbins)
			nbins = max(1, ceil((xmax - left_edge) / binwidth))
		else
			if nbins > 1
				ll = (xmax - left_edge) / nbins
				ul = (xmax - left_edge) / (nbins - 1)
				p10 = get_pow10(ul - ll)
				binwidth = p10 * ceil(ll / p10)
			end
		end
		right_edge = min(max(left_edge + nbins * binwidth, xmax), Float64(Inf))

	else
		error("TODO")
		# nbins = isnothing(nbins) ? 1 : nbins

		# binRange = max(1, ceil(nbins * eps(xscale)))
		# left_edge = floor(2 * (xmin - binRange / 4)) / 2
		# right_edge = ceil(2 * (xmax + binRange / 4)) / 2

		# binwidth = (right_edge - left_edge) / nbins
	end

	range(left_edge, stop = right_edge, length = nbins + 1)
end

function autorule(x::AbstractVector{Float64}, xmin::Float64, xmax::Float64, hardlimits::Bool)
	# Scott's normal reference rule
	binwidth = 3.5 * std(x) / (length(x)^(1 / 3))
	hardlimits ? error("TODO") : binpicker(xmin, xmax, nothing, binwidth)
end

function get_histcounts(
	x::AbstractVector{Float64};
	nbins::Union{Int64, Nothing} = nothing,
	binwidth::Union{Int64, Float64, Nothing} = nothing,
	norm::Symbol = :none,
	allow_nan = :false,
)
	if !allow_nan
		x = filter(x -> !isnan(x), x)
	end

	xmin, xmax = extrema(x)

	if isnothing(binwidth)
		binwidth = (xmax - xmin) / nbins
	end

	if isnothing(nbins)
		edges = autorule(x, xmin, xmax, false)
	else
		edges = binpicker(xmin, xmax, nbins, binwidth)
	end

	n, _ = histcountindices(x, edges)

	# normalization
	if norm != :none
		normalizations = Dict(
			:countdensity => () -> n ./ diff(edges, dims = 1),
			:cumcount => () -> cumsum(n),
			:probability => () -> n / length(x),
			:percentage => () -> (100 * n) / length(x),
			:pdf => () -> n / length(x) ./ diff(edges, dims = 1),
			:cdf => () -> cumsum(n / length(x)),
		)
		n = normalizations[norm]()
	end

	return n, edges
end

get_histcounts(x::AbstractVector{<:AbstractFloat}; kwargs...) = get_histcounts(Float64.(x); kwargs...)


