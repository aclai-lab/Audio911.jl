"""
	Audio911 signal data structures

	uses package Parameters for @with_kw mutable struct

	AudioSetup stores all datas that has to be shared in Audio911 module
	AudioData stores all results from signal analysis
"""
# ---------------------------------------------------------------------------- #
#                                    audio                                     #
# ---------------------------------------------------------------------------- #
struct Audio 
	data::AbstractVector{Float64}
	sr::Int64
end

# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
struct Stft
	# setup
	stft_length::Int64
	win::AbstractVector{Float64}
	win_type::Tuple{Symbol, Symbol}
	win_length::Int64
	overlap_length::Int64
	norm::Symbol # :none, :power, :magnitude, :pow2mag
	# data
	frames::AbstractArray{Float64}
	stft::AbstractArray{Float64}
	freq::AbstractVector{Float64}
end

# ---------------------------------------------------------------------------- #
#                             linear spectrogram                               #
# ---------------------------------------------------------------------------- #
struct LinSpec
	# setp
	win_norm::Symbol # :none, :power, :magnitude
	db_scale::Bool
	# data
	spec::AbstractArray{Float64}
	freq::AbstractVector{Float64}
end

# ---------------------------------------------------------------------------- #
#                                 filterbank                                   #
# ---------------------------------------------------------------------------- #
struct Fbank
	# setup
    bands::Int64
    scale::Symbol # :mel_htk, :mel_slaney, :erb, :bark
    norm::Symbol # :bandwidth, :area, :none
	# data
	fbank::AbstractArray{Float64}
	freq::AbstractVector{Float64}
end

# ---------------------------------------------------------------------------- #
#                                 audio rack                                   #
# ---------------------------------------------------------------------------- #

mutable struct AudioRack
	audio::Audio
	stft::Union{Stft, Nothing}
	lin::Union{LinSpec, Nothing}
	fbank::Union{Fbank, Nothing}

	# custom constructor
    function AudioRack(
        audio::AbstractVector{T},
        sr::Int64; 
        stft::Union{Stft, Nothing} = nothing,
        lin::Union{LinSpec, Nothing} = nothing,
        fb::Union{Fbank, Nothing} = nothing
    ) where T <: Real
        new(Audio(Float64.(audio), sr), stft, lin, fb)
    end
end

# @with_kw mutable struct AudioObj

# 	stft::Stft = Stft()
# 	lin::LinSpec = LinSpec()
# 	fb::Fbank = Fbank()

# 	# chroma
# 	bins_octave::Int64 # shared with constant-q transform
# 	center_freq::Int64
# 	gaussian_sd::Int64

# 	# mfcc
# 	mfcc_coeffs::Int64
# 	normalization_type::Symbol
# 	rectification::Symbol
# 	log_energy_source::Symbol
# 	log_energy_pos::Symbol
# 	delta_win_length::Int64
# 	delta_matrix::Symbol

# 	# spectral
# 	# spectrals::Spectrals
# 	spectral_spectrum::Symbol

# 	# f0
# 	f0_method::Symbol
# 	f0_range::Tuple{Int64, Int64}
# 	median_filter_length::Int64

# 	# constant-q transform
# 	freq_limits::Tuple{Float64, Float64}
# 	transform_type::Symbol
# end

# @with_kw mutable struct AudioData
# 	x::AbstractVector{Float64}

# 	stft::StftData = StftData()
# 	lin::LinData = LinData()
# 	fb::FbData = FbData()

# 	# mel_spectrum
# 	mel_frequencies::AbstractVector{Float64} = []
# 	mel_spectrogram::AbstractArray{Float64} = []

# 	# logaritmic mel
# 	log_mel::AbstractArray{Float64} = []

# 	# mfcc
# 	mfcc_coeffs::AbstractArray{Float64} = []
# 	mfcc_delta::AbstractArray{Float64} = []
# 	mfcc_deltadelta::AbstractArray{Float64} = []
# 	log_energy::AbstractVector{Float64} = []

# 	# spectral
# 	spectral_centroid::AbstractVector{Float64} = []
# 	spectral_crest::AbstractVector{Float64} = []
# 	spectral_decrease::AbstractVector{Float64} = []
# 	spectral_entropy::AbstractVector{Float64} = []
# 	spectral_flatness::AbstractVector{Float64} = []
# 	spectral_flux::AbstractVector{Float64} = []
# 	spectral_kurtosis::AbstractVector{Float64} = []
# 	spectral_rolloff::AbstractVector{Float64} = []
# 	spectral_skewness::AbstractVector{Float64} = []
# 	spectral_slope::AbstractVector{Float64} = []
# 	spectral_spread::AbstractVector{Float64} = []

# 	# f0
# 	f0::AbstractVector{Float64} = []

# 	# constant-q transform
# 	cqt_spec::AbstractArray{Float64} = []
# end