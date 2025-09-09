# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
abstract type AbstractStft end

# ---------------------------------------------------------------------------- #
#                                stft struct                                   #
# ---------------------------------------------------------------------------- #
struct Stft{F,T} <: AbstractStft
	stft_spec :: Matrix{T}
	stft_freq :: Vector{Float64}
	info      :: NamedTuple

	function Stft{F}(
		stft_spec :: Matrix{T},
		stft_freq :: Vector{Float64},
		info      :: NamedTuple
	) where {F,T}
		new{F,T}(stft_spec, stft_freq, info)
	end
end

#------------------------------------------------------------------------------#
#                                    methods                                   #
#------------------------------------------------------------------------------#
Base.eltype(::Stft{T}) where T = T

"""
    get_stft(s::Stft) -> Matrix{T}

Extract the STFT spectrogram matrix from an `Stft` container.

# See Also
- [`get_stft_freq`](@ref): Get frequency vector for the spectrogram
- [`get_info`](@ref): Get STFT computation parameters
- [`Stft`](@ref): STFT container type
"""
get_stft(s::Stft) = s.stft_spec

"""
    get_stft_freq(s::Stft) -> Vector{Float64}

Extract the frequency vector from an `Stft` container.

# See Also
- [`get_stft`](@ref): Get the STFT spectrogram matrix
- [`get_info`](@ref): Get STFT computation parameters  
- [`Stft`](@ref): STFT container type
"""
get_stft_freq(s::Stft) = s.stft_freq

"""
    get_info(s::Stft) -> NamedTuple

Extract  metadata from an `Stft` container.

# See Also
- [`get_stft`](@ref): Get the STFT spectrogram matrix
- [`get_stft_freq`](@ref): Get frequency vector
- [`Stft`](@ref): STFT container type
"""
get_info(s::Stft) = s.info

#------------------------------------------------------------------------------#
#                                  utilities                                   #
#------------------------------------------------------------------------------#
function get_onesided_fft_range(stft_size::Int64)
	if mod(stft_size, 2) == 0
		return collect(1:Int(stft_size / 2 + 1))   # EVEN
	else
		return collect(1:Int((stft_size + 1) / 2))  # ODD
	end
end # get_onesided_fft_range

#------------------------------------------------------------------------------#
#                                   get stft                                   #
#------------------------------------------------------------------------------#
function get_stft(
	frames          :: AudioFrames;
	stft_size       :: Int64=get_wsize(frames),
	frequency_range :: FreqRange=FreqRange(0, get_info(frames).sr÷2),
	spectrum_type   :: Symbol=:power, # :power, :magnitude
)::Stft
	sr = get_info(frames).sr
	win_size   = get_wsize(frames)
	overlap    = get_ovrlap(frames)
	n_channels = nchannels(frames)
	f          = reduce(hcat, get_frames(frames))

    # validate frequency range
    (0 ≤ frequency_range.low < frequency_range.hi ≤ sr / 2) ||
        throw(ArgumentError("Frequency range must be (0, sr÷2). " *
				"Got range: $(frequency_range.low), $(frequency_range.hi), sr÷2 = $(sr÷2)"))
    
    # validate overlap
    (0 < overlap < win_size) ||
        throw(ArgumentError("Overlap length must be < window length. " *
				"Got overlap = $overlap, window length = $win_size"))
    
    # validate stft_size
    stft_size < win_size &&
        throw(ArgumentError("stft_size must be ≥ window length. " *
				"Got stft_size = $stft_size, window length = $win_size"))
    
    # Validate mono audio
    n_channels == 1 ||
        throw(ArgumentError("STFT is designed to work on mono audiofiles only. " *
				"Got $n_channels channels"))

	# ensure frames is of length stft_size
	# if the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
	# this zero-padding in the time domain results in an interpolation in the frequency domain, 
	# which can provide a more detailed view of the spectral content of the signal.
	@inline @views f = win_size < stft_size ? 
		vcat(f, zeros(eltype(frames), stft_size - win_size, size(f, 2))) : 
		f[1:stft_size, :]

	# get fft
	stft_spec = fft(f, (1,))

	# post process
	# trim to desired range
	bin_low  = ceil(Int,  frequency_range.low * stft_size / sr + 1)
	bin_high = floor(Int, frequency_range.hi  * stft_size / sr + 1)
	bins     = collect(bin_low:bin_high)
	stft_spec = @views stft_spec[bins, :]

	# convert to half-sided power or magnitude spectrum
	spectrum_funcs = Dict(
		:power     => frames -> real.(frames .* conj.(frames)),
		:magnitude => frames -> abs.(frames),
	)
	# check if spectrum_type is valid
	@assert haskey(spectrum_funcs, spectrum_type) "Unknown spectrum_type: $spectrum_type."

	stft_spec = spectrum_funcs[spectrum_type](stft_spec)

	# trim borders
	# halve the first bin if it's the lowest bin
	bin_low  == 1 && (@views stft_spec[1, :] *= 0.5)
	# halve the last bin if it's the Nyquist bin and FFT length is even
	bin_high == fld(stft_size, 2) + 1 && iseven(stft_size) && (@views stft_spec[end, :] *= 0.5)

	# create frequency vector
	stft_freq = (sr / stft_size) * (bins .- 1)
	# shift final bin if fftLength is odd and the final range is full to fs/2.
	if stft_size % 2 != 0 && bin_high == floor(fftLength / 2 + 1)
		stft_freq[end] = sr * (stft_size - 1) / (2 * stft_size)
	end

	info = (;
		sr,
		stft_size,
		win_size,
		overlap,
		frequency_range,
		spectrum_type,
	)
	info = merge(info, get_info(frames))

	return Stft{AudioFrames}(stft_spec, stft_freq, info)
end

function get_stft(
	afile :: AudioFile;
    win   :: WinFunction=MovingWindow(
                            window_size=sr(afile)≤8000 ? 256 : 512,
                            window_step=sr(afile)≤8000 ? 128 : 256
                        ),
	type  :: Tuple{Symbol, Symbol}=(:hann, :periodic),
	kwargs...
	)::Stft
	frames = get_frames(afile; win, type)

	get_stft(frames; kwargs...)
end



# #------------------------------------------------------------------------------#
# #              fft version 1 as used in audio features extraction              #
# #------------------------------------------------------------------------------#
# function _get_fft(x::AbstractArray{Float64}, setup::AudioSetup)
# 	hop_length = setup.window_length - setup.get_ovrlap(frames)
# 	if isempty(setup.window)
# 		setup.window, _ = gencoswin(setup.window_type[1], setup.window_length, setup.window_type[2])
# 	end

# 	# split in windows
# 	y = buffer(x, setup.window_length, setup.window_length - setup.get_ovrlap(frames))

# 	# apply window and take fft
# 	Z = fft(y .* setup.window, (1,))

# 	# take one side
# 	logical_ossb = falses(setup.fft_length)
# 	logical_ossb[get_onesided_fft_range(setup.fft_length)] .= true
# 	Z = Z[logical_ossb, :]

# 	# log energy
# 	# reference: ETSI ES 201 108 V1.1.2 (2000-04)
# 	# https://www.3gpp.org/ftp/tsg_sa/TSG_SA/TSGS_13/Docs/PDF/SP-010566.pdf
# 	if setup.log_energy_pos != :none && setup.log_energy_source == :standard
# 		log_energy = sum(eachrow(y .^ 2))

# 		if setup.normalization_type == :standard
# 			log_energy[log_energy.==0] .= floatmin(Float64)
# 		elseif setup.normalization_type == :dithered
# 			log_energy[log_energy.<1e-8] .= 1e-8
# 		end
# 		log_energy = log.(log_energy)
# 	else
# 		log_energy = Vector{Float64}()
# 	end

# 	if setup.spectrum_type == :power
# 		real(Z .* conj(Z)), log_energy
# 	elseif setup.spectrum_type == :magnitude
# 		abs.(Z), log_energy
# 	else
# 		error("Unknown spectrum type: $(setup.spectrum_type)")
# 	end
# end

# get_fft!(setup::AudioSetup, data::AudioData) = data.fft, data.log_energy = _get_fft(data.x, setup)

# #------------------------------------------------------------------------------#
# #                                     stft                                     #
# #------------------------------------------------------------------------------#
# function _get_stft(
# 	x::AbstractArray{Float64},
# 	sr::Int64;
# 	fft_length::Int64,
# 	window::AbstractVector{Float64},
# 	frequency_range::Tuple{Int64, Int64},
# 	spectrum_type::Symbol,
# )
# 	@assert fft_length >= size(x, 1) "fft_length must be > window length. Got fft_length = $fft_length, window length = $(size(x,1))."

# 	# ensure x is of length fft_length
# 	# if the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
# 	# this zero-padding in the time domain results in an interpolation in the frequency domain, 
# 	# which can provide a more detailed view of the spectral content of the signal.
# 	x = size(x, 1) < fft_length ? vcat(x, zeros(eltype(x), fft_length - size(x, 1), size(x, 2))) : x[1:fft_length, :]

# 	# get fft
# 	Y = fft(x, (1,))

# 	# post process
# 	# trim to desired range
# 	bin_low = ceil(Int, frequency_range[1] * fft_length / sr + 1)
# 	bin_high = floor(Int, frequency_range[2] * fft_length / sr + 1)
# 	bins = collect(bin_low:bin_high)
# 	y = Y[bins, :]

# 	# convert to half-sided power or magnitude spectrum
# 	spectrum_funcs = Dict(
# 		:power => x -> real.((x .* conj.(x)) / (0.5 * sum(window.^2))),
# 		:magnitude => x -> abs.(x) / (0.5 * sum(window)),
# 	)
# 	# check if spectrum_type is valid
# 	@assert haskey(spectrum_funcs, spectrum_type) "Unknown spectrum_type: $spectrum_type."

# 	y = spectrum_funcs[spectrum_type](y)

# 	# trim borders
# 	# halve the first bin if it's the lowest bin
# 	bin_low == 1 && (y[1, :] *= 0.5)
# 	# halve the last bin if it's the Nyquist bin and FFT length is even
# 	bin_high == fld(fft_length, 2) + 1 && iseven(fft_length) && (y[end, :] *= 0.5)

# 	# create frequency vector
# 	stft_freq = (sr / fft_length) * (bins .- 1)
# 	# shift final bin if fftLength is odd and the final range is full to fs/2.
# 	if fft_length % 2 != 0 && bin_high == floor(fftLength / 2 + 1)
# 		stft_freq[end] = sr * (fft_length - 1) / (2 * fft_length)
# 	end

# 	return y, stft_freq
# end

# _get_stft(x::AbstractArray{Float64}, s::AudioSetup) = _get_stft(
# 	x,
# 	s.sr,
# 	fft_length = s.fft_length,
# 	window = s.window,
# 	frequency_range = s.frequency_range,
# 	spectrum_type = s.spectrum_type,
# )

# function get_stft(
# 	x::AbstractArray{<:AbstractFloat},
# 	sr::Int64,
# 	window::AbstractVector{Float64}=Float64[];
# 	fft_length::Int64 = sr <= 8000 ? 256 : 512,
# 	frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
# 	spectrum_type::Symbol = :power, # :power, :magnitude
# )
# 	@assert !isempty(window) "Window missing as third parameter."
# 	@assert sr > 0 "Sample rate must be > 0."
# 	@assert 0 <= frequency_range[1] < frequency_range[2] <= sr / 2 "Frequency range must be (0, sr/2)."

# 	stft_spec, stft_freq = _get_stft(
# 		x,
# 		sr,
# 		fft_length = fft_length,
# 		window = window,
# 		frequency_range = frequency_range,
# 		spectrum_type = spectrum_type,
# 	)
# end

# function get_stft(
# 	x::AbstractVector{<:AbstractFloat},
# 	sr::Int64;
# 	fft_length::Int64 = sr <= 8000 ? 256 : 512,
# 	window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
# 	window_length::Int64 = fft_length,
# 	get_ovrlap(frames)::Int64 = round(Int, fft_length / 2),
# 	frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
# 	spectrum_type::Symbol = :power, # :power, :magnitude
# )
# 	@assert sr > 0 "Sample rate must be > 0."
# 	@assert 0 <= frequency_range[1] < frequency_range[2] <= sr / 2 "Frequency range must be (0, sr/2)."
# 	@assert 0 < get_ovrlap(frames) < window_length "Overlap length must be < window length."

# 	frames, window = _get_frames(
# 		eltype(x) == Float64 ? x : Float64.(x),
# 		window_type = window_type,
# 		window_length = window_length,
# 		get_ovrlap(frames) = get_ovrlap(frames),
# 	)

# 	stft_spec, stft_freq = _get_stft(
# 		frames .* window,
# 		sr,
# 		fft_length = fft_length,
# 		window = window,
# 		frequency_range = frequency_range,
# 		spectrum_type = spectrum_type,
# 	)
# end

# function get_stft!(a::AudioObj)
# 	if isempty(a.data.frames)
# 		a.data.frames, a.setup.window = _get_frames(a.data.x, a.setup)
# 	end
# 	 a.data.stft, a.setup.stft_freq = _get_stft(a.data.frames .* a.setup.window, a.setup)
# end

# _get_stft2(x::AbstractArray{Float64}, s::AudioSetup) = _get_stft(
# 	x,
# 	s.sr,
# 	stft_length = s.stft_length,
# 	frequency_range = s.frequency_range,
# 	spectrum_type = s.spectrum_type,
# )

# function get_stft(
# 	frames::AudioFrames,
# 	sr::Int64;
# 	frequency_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
# 	spectrum_type::Symbol = :power, # :power, :magnitude
# )
# 	# @assert sr > 0 "Sample rate must be > 0."
# 	@assert 0 <= frequency_range[1] < frequency_range[2] <= sr / 2 "Frequency range must be (0, sr/2)."
# 	@assert 0 < get_ovrlap(frames) < window_length "Overlap length must be < window length."

# 	frames = _get_frames2(
# 		eltype(x) == Float64 ? x : Float64.(x),
# 		window_type = window_type,
# 		window_length = window_length,
# 		get_ovrlap(frames) = get_ovrlap(frames),
# 	)

# 	stft_spec, stft_freq = _get_stft2(
# 		frames .* window,
# 		sr,
# 		fft_length = fft_length,
# 		frequency_range = frequency_range,
# 		spectrum_type = spectrum_type,
# 	)
# end

