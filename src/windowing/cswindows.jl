function winhann(x)
	if x % 2 == 0
		xi = vcat(0:1/x:0.5-1/x, -0.5:1/x:-1/x)
	else
		xi = vcat(0:1/x:0.5-0.5/x, -0.5+0.5/x:1/x:-1/x)
	end
	g = 0.5 .+ 0.5 * cos.(2 * pi * xi)
	g = g .* (abs.(xi) .< 0.5)
	return g
end

function winhamm(x)
	if x % 2 == 0
		xi = vcat(0:1/x:0.5-1/x, -0.5:1/x:-1/x)
	else
		xi = vcat(0:1/x:0.5-0.5/x, -0.5+0.5/x:1/x:-1/x)
	end
	g = 0.54 .+ 0.46 * cos.(2 * pi * xi)
	g = g .* (abs.(xi) .< 0.5)
	return g
end

function winsine(x)
	if x % 2 == 0
		xi = vcat(0:1/x:0.5-1/x, -0.5:1/x:-1/x)
	else
		xi = vcat(0:1/x:0.5-0.5/x, -0.5+0.5/x:1/x:-1/x)
	end
	g = sin.(pi / 2 * cos.(pi * xi) .^ 2)
	g = g .* (abs.(xi) .< 0.5)
	return g
end

function winbh(x)
	if x % 2 == 0
		xi = vcat(0:1/x:0.5-1/x, -0.5:1/x:-1/x)
	else
		xi = vcat(0:1/x:0.5-0.5/x, -0.5+0.5/x:1/x:-1/x)
	end
	g = 0.35875 .+ 0.48829 * cos.(2 * pi * xi) .+ 0.14128 * cos.(4 * pi * xi) .+ 0.01168 * cos.(6 * pi * xi)
	g = g .* (abs.(xi) .< 0.5)
	return g
end

function wintriang(x)
	if x % 2 == 0
		xi = vcat(0:1/x:0.5-1/x, -0.5:1/x:-1/x)
	else
		xi = vcat(0:1/x:0.5-0.5/x, -0.5+0.5/x:1/x:-1/x)
	end
	g = 1 .- 2 * abs.(xi)
	g = g .* (abs.(xi) .< 0.5)
	return g
end

function cswindow(
	win_type::Symbol,
	bw::AbstractVector{Int64},
	nbins::Int64,
)
	nwin = length(bw)
	g = Array{Any}(undef, nwin)

	# winfun = getfield(Main, win_type)

	winfun = Dict([
		:hann => winhann,
		:hamming => winhamm,
		:itersine => winsine,
		:blackmanharris => winbh,
		:triang => wintriang])

	g = [winfun[win_type](bw[i]) for i in 1:nwin]

	if bw[1] > bw[2]
		g[1] = ones(bw[1], 1)
		g[1][(floor(Int, bw[1] / 2)-floor(Int, bw[2] / 2)+1):(floor(Int, bw[1] / 2)+ceil(Int, bw[2] / 2))] = winfun(bw[2])
	end

	if bw[nbins+2] > bw[nbins+3]
		g[nbins+2] = ones(bw[nbins+2], 1)
		g[nbins+2][(floor(Int, bw[nbins+2] / 2)-floor(Int, bw[nbins+3] / 2)+1):(floor(Int, bw[nbins+2] / 2)+ceil(Int, bw[nbins+3] / 2))] = winfun(bw[nbins+3])
	end

    return g
end

