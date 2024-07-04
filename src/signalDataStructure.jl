"""
	Audio911 signal data structures

	uses package Parameters for @with_kw mutable struct

	AudioSetup stores all datas that has to be shared in Audio911 module
	AudioData stores all results from signal analysis
"""
# ---------------------------------------------------------------------------- #
#                                 raw audio                                    #
# ---------------------------------------------------------------------------- #
@with_kw mutable struct Audio
	audio::AbstractVector{Float64} = Float64[]
	sr::Int64
end

# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
@with_kw mutable struct Stft
	# setup
	stft_length::Int64
	win::AbstractVector{Float64} = Float64[]
	win_type::Tuple{Symbol, Symbol}
	win_length::Int64
	overlap_length::Int64
	norm::Symbol # :none, :power, :magnitude, :pow2mag
	# data
	frames::AbstractArray{Float64} = Float64[]
	stft::AbstractArray{Float64} = Float64[]
	freq::AbstractVector{Float64} = Float64[]
end

# ---------------------------------------------------------------------------- #
#                             linear spectrogram                               #
# ---------------------------------------------------------------------------- #
@with_kw mutable struct LinSpec
	# setp
	win_norm::Symbol # :none, :power, :magnitude
	db_scale::Bool
	# data
	spec::AbstractArray{Float64} = Float64[]
	freq::AbstractVector{Float64} = Float64[]
end

# ---------------------------------------------------------------------------- #
#                                 filterbank                                   #
# ---------------------------------------------------------------------------- #
@with_kw mutable struct Fbank
	# setup
    bands::Int64
    scale::Symbol # :mel, :erb, :bark
    norm::Symbol # :bandwidth, :area, :none
    mel_style::Symbol # :htk, :slaney
	# data
	fbank::AbstractArray{Float64} = []
	freq::AbstractVector{Float64} = []
end

# ---------------------------------------------------------------------------- #
#                                audio object                                  #
# ---------------------------------------------------------------------------- #
@with_kw mutable struct AudioObj
	audio::Audio
	stft::Stft
	lin::LinSpec
	fb::Fbank

	# chroma
	bins_octave::Int64 # shared with constant-q transform
	center_freq::Int64
	gaussian_sd::Int64

	# mfcc
	mfcc_coeffs::Int64
	normalization_type::Symbol
	rectification::Symbol
	log_energy_source::Symbol
	log_energy_pos::Symbol
	delta_win_length::Int64
	delta_matrix::Symbol

	# spectral
	# spectrals::Spectrals
	spectral_spectrum::Symbol

	# f0
	f0_method::Symbol
	f0_range::Tuple{Int64, Int64}
	median_filter_length::Int64

	# constant-q transform
	freq_limits::Tuple{Float64, Float64}
	transform_type::Symbol
end

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