"""
	Audio911 signal data structures

	uses package Parameters for @with_kw mutable struct

	AudioSetup stores all datas that has to be shared in Audio911 module
	AudioData stores all results from signal analysis
"""
@with_kw mutable struct StftSetup
	stft_length::Int64
	# windowing
	win::AbstractVector{Float64} = []
	win_type::Tuple{Symbol, Symbol}
	win_length::Int64
	overlap_length::Int64
	# spectrum
	spec_norm::Symbol # :none, :power, :magnitude, :pow2mag
end

@with_kw mutable struct StftData
	frames::AbstractArray{Float64} = []
	stft::AbstractArray{Float64} = []
	freq::AbstractVector{Float64} = []
end

@with_kw mutable struct LinSetup
	win_norm::Symbol # :none, :power, :magnitude
	db_scale::Bool
end

@with_kw mutable struct LinData
	spec::AbstractArray{Float64} = []
	freq::AbstractVector{Float64} = []
end

# @with_kw mutable struct FbSetup
# 	n_bands::Int64
# 	design_domain::Symbol, # :linear, :mel, :erb, :bark
# 	mel_style::Symbol
# end

# @with_kw mutable struct FbData

# end

@with_kw mutable struct AudioSetup
	sr::Int64
	freq_range::Tuple{Int64, Int64}

	stft::StftSetup
	lin::LinSetup

	# mel
	
	mel_style::Symbol
	mel_bands::Int64
	design_domain::Symbol
	fb_norm::Symbol
	frequency_scale::Symbol
	st_peak_range::Tuple{Int64, Int64}

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

@with_kw mutable struct AudioData
	x::AbstractVector{Float64}

	stft::StftData = StftData()
	lin::LinData = LinData()

	# mel_spectrum
	mel_frequencies::AbstractVector{Float64} = []
	mel_spectrogram::AbstractArray{Float64} = []

	# logaritmic mel
	log_mel::AbstractArray{Float64} = []

	# mfcc
	mfcc_coeffs::AbstractArray{Float64} = []
	mfcc_delta::AbstractArray{Float64} = []
	mfcc_deltadelta::AbstractArray{Float64} = []
	log_energy::AbstractVector{Float64} = []

	# spectral
	spectral_centroid::AbstractVector{Float64} = []
	spectral_crest::AbstractVector{Float64} = []
	spectral_decrease::AbstractVector{Float64} = []
	spectral_entropy::AbstractVector{Float64} = []
	spectral_flatness::AbstractVector{Float64} = []
	spectral_flux::AbstractVector{Float64} = []
	spectral_kurtosis::AbstractVector{Float64} = []
	spectral_rolloff::AbstractVector{Float64} = []
	spectral_skewness::AbstractVector{Float64} = []
	spectral_slope::AbstractVector{Float64} = []
	spectral_spread::AbstractVector{Float64} = []

	# f0
	f0::AbstractVector{Float64} = []

	# constant-q transform
	cqt_spec::AbstractArray{Float64} = []
end

mutable struct AudioObj
	setup::AudioSetup
	data::AudioData
end