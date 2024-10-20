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
dotu(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = dot(x, y)
function dotu(x::AbstractVector{T}, y::AbstractVector{V}) where {T,V}
    dotprod = zero(promote_type(T, V))
    for i in eachindex(x, y)
        dotprod = muladd(x[i], y[i], dotprod)
    end
    dotprod
end

function levinson(R_xx::AbstractVector{U}, p::Integer) where U<:Number
    # for m = 1
    k = -R_xx[2] / R_xx[1]
    F = promote_type(Base.promote_union(U), typeof(k))
    prediction_err = real(R_xx[1] * (one(F) - abs2(k)))
    R = typeof(prediction_err)
    T = promote_type(F, R)

    a = zeros(T, p)
    reflection_coeffs = zeros(T, p)
    a[1] = reflection_coeffs[1] = k
    rev_a = similar(a, p - 1)     # buffer to store a in reverse

    @views for m = 2:p
        rev_a[1:m-1] .= a[m-1:-1:1]
        k = -(R_xx[m+1] + dotu(R_xx[2:m], rev_a[1:m-1])) / prediction_err
        @. a[1:m-1] = muladd(k, conj(rev_a[1:m-1]), a[1:m-1])
        a[m] = reflection_coeffs[m] = k
        prediction_err *= (one(R) - abs2(k))
    end

    # Return autocorrelation coefficients, error estimate, and reflection coefficients
    return a, prediction_err, reflection_coeffs
end

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

function pitch_estimation_filter(freq::AbstractVector{<:AbstractFloat})
	k = 10
	gamma = 1.8
	num = round(Int, length(freq) / 2)
	q = 10 .^ range(log10(0.5), log10(k + 0.5), length=num)
	h = @. 1 / (gamma - cos(2π * q))
	delta = diff([q[1]; (q[1:end-1] .+ q[2:end]) ./ 2; q[end]])
	beta = sum(h .* delta) / sum(delta)

	a_filt = (h .- beta)'
	npad = findfirst(x -> x >= 1, q) - 1

	return a_filt, npad
end

function get_trimmedspec(x::AbstractArray{<:AbstractFloat}, nfft::Int64, maxbin::Int64)
	win = _get_window(:hamming, size(x, 1), true)
	x_reframed = _get_wframes(x, win)
	x_pad = (vcat(col, zeros(eltype(col), nfft - nwin)) for col in x_reframed)
	fft(combinedims(collect(x_pad)), 1)[1:maxbin, :]
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
		win = _get_window(:blackman, nwin, true)
		residual_buff = _get_wframes(combinedims(residual_buff), win)
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

	elseif method == :pef
		nfft = nextpow(2, 2 * nwin - 1)
		ncol = size(x, 2)
		logspaced_freq = 10 .^ range(1, log10(min(sr / 2 - 1, 4000)), length=nfft)
		linspaced_freq = range(0, sr / 2, length=round(Int, nfft / 2) + 1)

		logfrqrange = map(x -> findmin(abs.(logspaced_freq .- x))[2], freqrange)

		bwtmp = (logspaced_freq[3:end] - logspaced_freq[1:end-2]) / 2
		bw = [bwtmp[1]; bwtmp; bwtmp[end]] ./ nfft
		
		a_filt, npad = pitch_estimation_filter(logspaced_freq)
		
		win = _get_window(:hamming, nwin, true)
		x_reframed = _get_wframes(x, win)
		x_pad = (vcat(col, zeros(eltype(col), nfft - nwin)) for col in x_reframed)
		x_refft = fft(combinedims(collect(x_pad)), 1)[1:div(nfft, 2)+1, :]
		x_refft = real.(x_refft .* conj(x_refft))
		x_log = hcat([LinearInterpolation(linspaced_freq, col)(logspaced_freq) for col in eachcol(x_refft)]...) .* bw
		z = [zeros(eltype(x_log), npad, size(x_log, 2)); x_log]

		m = max(size(z, 1), size(a_filt, 1))
		mxl = min(logfrqrange[2], m - 1)
		m2 = min(nextpow(2, 2*m - 1), nfft * 4)

		z_pad = combinedims(collect(vcat(col, zeros(eltype(col), m2 - size(z, 1))) for col in eachcol(z)))
		a_filt_pad = vcat(a_filt', zeros(eltype(a_filt), m2 - length(a_filt)))
		c1 = real.(ifft(fft(z_pad, 1) .* conj.(fft(a_filt_pad)), 1))

		r = [c1[m2 .- mxl .+ (1:mxl), :]; c1[1:mxl+1, :]]

		domain = r[logfrqrange[2]+1:end, :]

		_, locs = get_candidates(domain, logfrqrange)

		f0 = logspaced_freq[vec(locs)]

	elseif method == :cep
		nfft = nextpow(2, 2 * nwin - 1)
		win = _get_window(:hamming, nwin, true)
		x_reframed = _get_wframes(x, win)
		x_pad = (vcat(col, zeros(eltype(col), nfft - nwin)) for col in x_reframed)
		x_refft = fft(combinedims(collect(x_pad)), 1)[1:div(nfft, 2)+1, :]

		domain = real.(ifft(log.(abs.(fft(combinedims(collect(x_pad)), 1)).^2), 1))

		_, locs = get_candidates(domain, freqrange)

		f0 = sr ./ vec(locs)

	elseif method == :lhs
		maxbin = 5 * (freqrange[2] + 1)
		s = log.(abs.(get_trimmedspec(x, sr, maxbin)))

		@views domain = 
			s[1:(freqrange[2]+1), :] + 
			s[1:2:(freqrange[2]+1)*2, :] +
			s[1:3:(freqrange[2]+1)*3, :] +
			s[1:4:(freqrange[2]+1)*4, :] +
			s[1:5:end, :]

		_, locs = get_candidates(domain, freqrange)
		f0 = vec(Float64.(locs))

	else
		error("Method `$method` not implemented.")
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
