# using FFTW
# using Audio911

# include("../windowing/cswindows.jl")

function cfs_coefs(cfs, n, sigType)
	# Multiplication of cfs by scaling factors
	isMATLAB = isempty(coder.target)
	if lowercase(sigType) == "real"

		cfs = [2 * size(x, 1) / n .* x for x in cfs]

	else

		cfs = [size(x, 1) / n .* x for x in cfs]

	end
end

function cqt_full(
	xDFT::AbstractVector{ComplexF64},
	g::AbstractVector{Vector{Float64}},
	cfbins::AbstractVector{Float64},
	m::AbstractVector{Float64},
	transform_type::Symbol,
	nyq_bin::Int64,
)

	nwin = length(g)
	ndft = length(xDFT)
	fIntervals = Vector{Vector{Int}}(undef, nwin)
	tCFS = Vector{Vector{ComplexF64}}(undef, nwin)

	# Algorithm for forward CQT due to Holighaus and Velasco
	for kk in 1:nwin
		g_length = length(g[kk])
		win_order = vcat([ceil(Int, g_length / 2)+1:g_length, 1:ceil(Int, g_length / 2)]...)
		win_range = Int.(1 .+ mod.(cfbins[kk] .+ (-floor(g_length / 2):ceil(g_length / 2)-1), ndft))
		fIntervals[kk] = win_range
		rowDim = Int(m[kk])

		tmp = zeros(ComplexF64, Int(m[kk]))
		tmp_indexes = vcat([rowDim-floor(Int, g_length / 2)+1:rowDim, 1:ceil(Int, g_length / 2)]...)
		tmp[tmp_indexes] = xDFT[win_range, :] .* g[kk][win_order]

		tCFS[kk] = ifft(tmp)
	end

	# function cfsCoefs
	tCFS = [2 * size(x, 1) / ndft .* x for x in tCFS]

	if transform_type == :reduced

		error("TODO!")
		# DCcfs = tCFS[1]
		# Nyquistcfs = tCFS[nyq_bin]

		# tCFS = tCFS[[2:nyq_bin-1; nyq_bin+1:end]]
		# tCFSMatrix = reduce(hcat, tCFS)

		# numTPts = maximum(m[2:nyq_bin-1])
		# tCFSReshape = reshape(tCFSMatrix, numTPts, nWin - 2, numSig)
		# # Permute frequency and hop so that the coefficients matrices are
		# # frequency by time
		# c = permutedims(tCFSReshape, [2, 1, 3])
		# cfsReduced = (c = c, DCcfs = DCcfs, Nyquistcfs = Nyquistcfs, NyquistBin = nyq_bin)
		# cfsFull = nothing
		return nothing
	else
		cfs_full = reduce(hcat, tCFS)'
		return cfs_full, fIntervals
	end
end

# win type :hann, :hamm, :itersine, :blackmanharris, :triang
# transform type :full, :reduced, :sparse
function cqt(
	x::AbstractVector{Float64},
	sr::Int64,
	nbins_octave::Int64 = 12,
	freq_limits::Tuple{Float64, Float64} = (0.0, 0.0),
	transform_type::Symbol = :full,
)
	x_length = length(x)

	if x_length < 4
		error("signal too small: less than 4 samples.")
	end

	if freq_limits == (0, 0)
		freq_limits = (sr / x_length, sr / 2)
	end

	# Nyquist frequency
	nyqfreq = sr / 2

	# Set minimum bandwidth in DFT bins
	minbw = 4

	q = 1 / (2^(1 / nbins_octave) - 2^(-(1 / nbins_octave)))
	# bw = cf*q^(-1) so obtain 1/q
	q_inv = 1 / q
	# Default number of octaves. Velasco et al., 2011 CQ-NSGT Parameters:
	# Windows and Lattices
	b_max = ceil(Int, nbins_octave * log2(freq_limits[2] / freq_limits[1]) + 1)

	tfreq_bins = freq_limits[1] .* 2 .^ ((0:b_max-1) ./ nbins_octave)

	id_max = findfirst(tfreq_bins .>= freq_limits[2])

	# Frequency bins less than Nyquist frequency
	if tfreq_bins[id_max] >= nyqfreq
		nyq_bins = tfreq_bins[1:id_max-1]
	else
		nyq_bins = tfreq_bins[1:id_max]
	end

	# First record number of bins, 1,....K prior to pre-pending DC and appending the Nyquist
	nbins = length(nyq_bins)
	# Now prepend DC and append Nyquist frequencies
	t_freq = vcat(0, nyq_bins, nyqfreq)
	# Store Nyquist Bin
	nyq_bin = nbins + 2
	# Store DC bin
	dc_bin = 1
	# Mirror other filters -- start with one bin below the Nyquist bin and go down to one bin above DC
	f = vcat(t_freq..., Float64(sr) .- t_freq[end-1:-1:2]...)

	# Convert the frequency bins to approximate index of DFT bin.
	fbins = f * (x_length / sr)

	# Determine bandwidths in DFT bins. For everything but the DC bin and the
	# Nyquist bins the bandwidth is \epsilon_{k+1}-\epsilon_{k-1}. Total number
	# of bandwidths is now 2*lenBins+2
	bw = zeros(2 * nbins + 2)

	#   Set bandwidth of DC bin to 2*fMin -- these are bandwidths in samples
	#   (approximately), we will round to integer values.
	bw[dc_bin] = 2 * fbins[2]
	# Set lenBins 1 such that cf/bw = q
	bw[2] = fbins[2] * q_inv
	#   Set the bandwidth for the frequency before the Nyquist such that
	#   cf/bw = q
	bw[nbins+1] = fbins[nbins+1] * q_inv
	# Set the original k = 1,....K-1
	idxk = vcat(3:nbins, nyq_bin)
	# See Velasco et al. and Holighaus et al. CQ-NSGT Parameters: Windows
	# and Lattices
	bw[idxk] = fbins[idxk.+1] - fbins[idxk.-1]
	# Mirror bandwidths on the negative frequencies
	bw[nbins+3:2*nbins+2] = bw[nbins+1:-1:2]
	# Round the bandwidths to integers
	bw = round.(Int, bw)

	# Convert frequency centers to integers. Round down up to Nyquist. Round up
	# after Nyquist.
	cfbins = zeros(size(fbins))
	# Up to Nyquist floor round down
	cfbins[1:nbins+2] = floor.(Int, fbins[1:nbins+2])
	# From Nyquist to Fs, round up
	cfbins[nbins+3:end] = ceil.(Int, fbins[nbins+3:end])

	# Compute the shift between filters in frequency in samples
	diffLFDC = x_length - cfbins[end]
	fshifts = vcat(diffLFDC, diff(cfbins))

	# Ensure that any bandwidth less than the minimum window is set to the minimum window
	bw[bw.<minbw] .= minbw

	# Compute the frequency windows for the CQT-NSGFT
	win_type = :hann # default window type TODO parametrization?
	g = cswindow(win_type, bw, nbins)

	# Obtain DFT of input signal
	xDFT = fft(x)

	# Depending on transform type value
	if transform_type == :full
		m = maximum(length.(g)) * ones(length(g))
		cfs, fintervals = cqt_full(xDFT, g, cfbins, m, transform_type, nyq_bin)

	elseif transform_type == :reduced
		error("TODO!")
		# m = zeros(size(bw))
		# hfBand = bw[nyq_bin-1]
		# m[2:nyq_bin, nyq_bin+1:end] = hfBand
		# m[dcBin] = bw[dcBin]
		# m[nyq_bin] = bw[nyq_bin]
		# _, cfs, fintervals = cqt_full(xDFT, g, cfbins, m, transform_type, nyq_bin)

	elseif transform_type == :sparse
		error("TODO!")
		# cfs, fintervals = cqt_sparse(xDFT, g, cfbins, bw)

	else
		error("Unrecognized string choice for transform_type")
	end

	return cfs, fintervals, nyq_bin, f
end

function get_cqt_spec!(
	setup::AudioSetup,
	data::AudioData,
)
	cqt_cspec, _, nyq_bin, f = cqt(data.x, setup.sr, setup.bins_octave, setup.freq_limits, setup.transform_type)

	if setup.transform_type == :full || setup.transform_type == :sparse
		f = f[1:nyq_bin]
		cqt_cspec = cqt_cspec[1:nyq_bin, :]
	elseif setup.transform_type == :reduced
		error("TODO!")
		# f = f[2:nyq_bin - 1]
		# # Nyquist frequency has already been removed from c field
		# cqt_cspec = cqt_cspec.c[1:nyq_bin - 2,:]
	else
		error("Unknown transform_type.")
	end

	data.cqt_spec = 20 * log10.(abs.(cqt_cspec) .+ eps(0.0))
end