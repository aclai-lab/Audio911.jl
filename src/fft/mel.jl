# ---------------------------------------------------------------------------- #
#                         mel scale utility functions                          #
# ---------------------------------------------------------------------------- #
function hz2mel(
	hz::Tuple{Int64, Int64},
	mel_style::Symbol = :htk, # :htk, :slaney
)
	if mel_style == :htk
		mel = 2595 * log10.(1 .+ reduce(vcat, getindex.(hz)) / 700)
	else # slaney
		hz = reduce(vcat, getindex.(hz))
		linStep = 200 / 3
		logStep = log(6.4) / 27
		changePoint = 1000
		changePoint_mel = changePoint / linStep
		isLinearRegion = hz .< changePoint
		mel = Float64.(hz)
		mel[isLinearRegion] .= hz[isLinearRegion] / linStep
		mel[.!isLinearRegion] .= changePoint_mel .+
								 log.(hz[.!isLinearRegion] / changePoint) / logStep
	end
	return mel
end # hz2mel

function mel2hz(
	mel::LinRange{Float64, Int64},
	mel_style::Symbol = :htk, # :htk, :slaney
)
	if mel_style == :htk
		hz = 700 * (exp10.(mel / 2595) .- 1)
	else
		linStep = 200 / 3
		logStep = log(6.4) / 27
		changePoint = 1000
		changePoint_mel = changePoint / linStep
		isLinearRegion = mel .< changePoint_mel
		hz = [mel;]
		hz[isLinearRegion] .= hz[isLinearRegion] * linStep
		hz[.!isLinearRegion] .= changePoint *
								exp.(logStep * (mel[.!isLinearRegion] .- changePoint_mel))
	end
	return hz
end # mel2hz

function get_mel_norm_factor(spectrum_type::Symbol, fft_window::Vector{Float64})
	if spectrum_type == :power
		return 1 / (sum(fft_window)^2)
	elseif spectrum_type == :magnitude
		return 1 / sum(fft_window)
	else
		error("Unknown spectrum_type $spectrum_type.")
	end
end

# ---------------------------------------------------------------------------- #
#                       semitones scale utility functions                      #
# ---------------------------------------------------------------------------- #
hz2semitone(freq, base_freq) = 12 * log2(freq / base_freq)
semitone2hz(z, base_freq)    = base_freq * (2^(z / 12))

function catch_loudest_index(
	lin_spec::AbstractArray{Float64},
	freq_array::AbstractVector{Float64},
	st_peak_range::Tuple{Int64, Int64},
)
	comp_spec = mean(lin_spec, dims = 1)'

	min_index = findfirst((x) -> x >= st_peak_range[1], freq_array)
	max_index = findfirst((x) -> x >= st_peak_range[2], freq_array) - 1

	_, loudest_index = findmax(comp_spec[min_index:max_index])

	return min_index + loudest_index
end
# ---------------------------------------------------------------------------- #
#                           design filterbank matrix                           #
# ---------------------------------------------------------------------------- #
### da generalizzare per 
### frequency_scale :mel, :bark, :erb
### filterbanl_design_domain :linear, :warped (da verificare se serve)
function design_filterbank(data::AudioData, setup::AudioSetup)
	# set the design domain ### da implementare in futuro
	setup.filterbank_design_domain == :linear ? design_domain = :linear :
	design_domain = setup.frequency_scale

	# compute band edges
	# TODO da inserire il caso :erb e :bark

	melRange = hz2mel(setup.frequency_range, setup.mel_style)

	# mimic audioflux linear mel_style
	if setup.mel_style == :linear
		lin_fq = collect(0:(setup.fft_length-1)) / setup.fft_length * setup.sr
		band_edges = lin_fq[1:(setup.mel_bands+2)]
	elseif setup.mel_style == :htk || setup.mel_style == :slaney
		band_edges = mel2hz(LinRange(melRange[1], melRange[end], setup.mel_bands + 2), setup.mel_style)

	elseif setup.mel_style == :semitones
		if isempty(data.lin_spectrogram)
			lin_spectrogram!(setup, data)
		end

		loudest_index = catch_loudest_index(data.lin_spectrogram, data.lin_frequencies, setup.st_peak_range)

		minsemitone = hz2semitone(setup.st_peak_range[1], loudest_index)
		maxsemitone = hz2semitone(setup.st_peak_range[2], loudest_index)

		band_edges = semitone2hz.(minsemitone .+ collect(0:(setup.mel_bands+1)) / (setup.mel_bands + 1) * (maxsemitone - minsemitone), loudest_index)

		println(band_edges)
	else
		error("Unknown mel_style $(setup.mel_style).")
	end

	### parte esclusiva per mel filterbank si passa a file designmelfilterbank.m
	# determine the number of bands
	num_edges = length(band_edges)

	# determine the number of valid bands
	valid_num_edges = sum((band_edges .- (setup.sr / 2)) .< sqrt(eps(Float64)))
	valid_num_bands = valid_num_edges - 2

	# preallocate the filter bank
	filterbank = zeros(Float64, setup.fft_length, setup.mel_bands)
	data.mel_frequencies = band_edges[2:(end-1)]

	# Set this flag to true if the number of FFT length is insufficient to
	# compute the specified number of mel bands
	FFTLengthTooSmall = false

	# if :hz 
	linFq = collect(0:(setup.fft_length-1)) / setup.fft_length * setup.sr

	# Determine inflection points
	@assert(valid_num_edges <= num_edges)
	p = zeros(Float64, valid_num_edges, 1)

	for edge_n in 1:valid_num_edges
		for index in eachindex(linFq)
			if linFq[index] > band_edges[edge_n]
				p[edge_n] = index
				break
			end
		end
	end

	FqMod = linFq

	# Create triangular filters for each band
	bw = diff(band_edges)

	for k in 1:Int(valid_num_bands)
		# Rising side of triangle
		for j in Int(p[k]):(Int(p[k+1])-1)
			filterbank[j, k] = (FqMod[j] - band_edges[k]) / bw[k]
		end
		# Falling side of triangle
		for j in Int(p[k+1]):(Int(p[k+2])-1)
			filterbank[j, k] = (band_edges[k+2] - FqMod[j]) / bw[k+1]
		end
		emptyRange1 = p[k] .> p[k+1] - 1
		emptyRange2 = p[k+1] .> p[k+2] - 1
		if (!FFTLengthTooSmall && (emptyRange1 || emptyRange2))
			FFTLengthTooSmall = true
		end
	end

	# mirror two sided
	range = get_onesided_fft_range(setup.fft_length)
	range = range[2:end]
	filterbank[end:-1:(end-length(range)+1), :] = filterbank[range, :]

	filterbank = filterbank'

	# normalizzazione    
	BW = band_edges[3:end] - band_edges[1:(end-2)]

	if (setup.filterbank_normalization == :area)
		weight_per_band = sum(filterbank, dims = 2)
		if setup.frequency_scale != :erb
			weight_per_band = weight_per_band / 2
		end
	elseif (setup.filterbank_normalization == :bandwidth)
		weight_per_band = BW / 2
	else
		weight_per_band = ones(1, setup.mel_bands)
	end

	for i in 1:(setup.mel_bands)
		if (weight_per_band[i] != 0)
			filterbank[i, :] = filterbank[i, :] ./ weight_per_band[i]
		end
	end

	# get one side
	range = get_onesided_fft_range(setup.fft_length)
	filterbank = filterbank[:, range]
	# manca la parte relativa a :erb e :bark

	# setta fattore di normalizzazione
	if setup.window_norm
		win_norm_factor = get_mel_norm_factor(setup.spectrum_type, data.fft_window)
		filterbank = filterbank * win_norm_factor
	end

	return filterbank
end # function design_filterbank

# ---------------------------------------------------------------------------- #
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
function get_mel_spec!(
	setup::AudioSetup,
	data::AudioData,
)
	filterbank = design_filterbank(data, setup)

	hop_length = setup.window_length - setup.overlap_length
	num_hops = Int(floor((size(data.x, 1) - setup.window_length) / hop_length) + 1)

	# apply filterbank
	# if (setup.spectrum_type == :power)
	data.mel_spectrogram = reshape(
		filterbank * data.fft, setup.mel_bands, num_hops)
	# else
	#     #TODO
	#     error("magnitude not yet implemented.")
	# end

	data.mel_spectrogram = transpose(data.mel_spectrogram)
end # melSpectrogram

# ---------------------------------------------------------------------------- #
#                          logaritmic mel spectrogram                          #
# ---------------------------------------------------------------------------- #
# TODO prova a fare le delta del log mel
function get_log_mel!(
	setup::AudioSetup,
	data::AudioData,
)
	# Reference:
	# https://dsp.stackexchange.com/questions/85501/log-of-filterbank-energies

	# Rectify
	mel_spec = deepcopy(data.mel_spectrogram')

	if setup.normalization_type == :standard
		mel_spec[mel_spec.==0] .= floatmin(Float64)
	elseif setup.normalization_type == :dithered
		mel_spec[mel_spec.<1e-8] .= 1e-8
	else
		@warn("Unknown $setup.normalization_type normalization type, defaulting to standard.")
		mel_spec[mel_spec.==0] .= floatmin(Float64)
	end

	data.log_mel = log10.(mel_spec)'
end

# ---------------------------------------------------------------------------- #
#                             mfcc related functions                           #
# ---------------------------------------------------------------------------- #
function create_DCT_matrix(
	mel_coeffs::Int64,
)
	# create DCT matrix
	matrix = zeros(Float64, mel_coeffs, mel_coeffs)
	s0 = sqrt(1 / mel_coeffs)
	s1 = sqrt(2 / mel_coeffs)
	piCCast = 2 * pi / (2 * mel_coeffs)

	matrix[1, :] .= s0
	for k in 1:mel_coeffs, n in 2:mel_coeffs
		matrix[n, k] = s1 * cos(piCCast * (n - 1) * (k - 0.5))
	end

	matrix
end

function audioDelta(
	x::AbstractMatrix{T},
	window_length::Int64,
	source::Symbol = :standard,
) where {T <: AbstractFloat}

	# define window shape
	m = Int(floor(window_length / 2))
	b = collect(m:-1:(-m)) ./ sum((1:m) .^ 2)

	if source == :transposed
		filt(b, 1.0, x')'   #:audioflux setting
	else
		filt(b, 1.0, x)     #:matlab setting
	end
end

function get_mfcc!(
	setup::AudioSetup,
	data::AudioData,
)
	# Rectify
	mel_spec = deepcopy(data.mel_spectrogram')

	if setup.normalization_type == :standard
		mel_spec[mel_spec.==0] .= floatmin(Float64)
	elseif setup.normalization_type == :dithered
		mel_spec[mel_spec.<1e-8] .= 1e-8
	else
		@warn("Unknown $setup.normalization_type normalization type, defaulting to standard.")
		mel_spec[mel_spec.==0] .= floatmin(Float64)
	end

	# Design DCT matrix
	DCTmatrix = create_DCT_matrix(setup.mel_bands)

	# apply DCT matrix
	if (setup.rectification == :log)
		coeffs = DCTmatrix * log10.(mel_spec)
	elseif (setup.rectification == :cubic_root)
		# apply DCT matrix
		coeffs = DCTmatrix * mel_spec .^ (1 / 3)
	else
		@warn("Unknown $rectification DCT matrix rectification, defaulting to log.")
		coeffs = DCTmatrix * log10.(mel_spec)
	end

	# reduce to mfcc coefficients
	data.mfcc_coeffs = coeffs[1:(setup.num_coeffs), :]'

	# log energy calc
	if setup.log_energy_source == :mfcc
		log_energy = sum(eachrow(mel_spec .^ 2)) / setup.mel_bands

		if setup.normalization_type == :standard
			log_energy[log_energy.==0] .= floatmin(Float64)
		elseif setup.normalization_type == :dithered
			log_energy[log_energy.<1e-8] .= 1e-8
		end

		data.log_energy = log.(log_energy)
	end

	if (setup.log_energy_pos == :append)
		data.mfcc_coeffs = hcat(data.mfcc_coeffs, data.log_energy)
	elseif (setup.log_energy_pos == :replace)
		data.mfcc_coeffs = hcat(data.log_energy, data.mfcc_coeffs[:, 2:end])
	end
end

function get_mfcc_deltas!(
	setup::AudioSetup,
	data::AudioData,
)
	data.mfcc_delta = audioDelta(
		data.mfcc_coeffs, setup.delta_window_length, setup.delta_matrix)
	data.mfcc_deltadelta = audioDelta(
		data.mfcc_delta, setup.delta_window_length, setup.delta_matrix)
end
