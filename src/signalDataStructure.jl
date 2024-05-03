"""
	Audio911 signal data structures

	uses package Parameters for @with_kw mutable struct

	AudioSetup stores all datas that has to be shared in Audio911 module
	AudioData stores all results from signal analysis
"""
@with_kw mutable struct AudioSetup
	sr::Int64

	# fft
	fft_length::Int64
	window_type::Tuple{Symbol, Symbol}
	window_length::Int64
	overlap_length::Int64
	window_norm::Bool

	# spectrum
	frequency_range::Tuple{Int64, Int64}
	lin_frequencies::Vector{Float64} = []
	band_edges::AbstractVector{AbstractFloat} = []
	spectrum_type::Symbol

	# mel
	mel_style::Symbol
	mel_bands::Int64
	mel_frequencies::Vector{Float64} = []
	filterbank_design_domain::Symbol
	filterbank_normalization::Symbol
	frequency_scale::Symbol

	# mfcc
	num_coeffs::Int64
	normalization_type::Symbol
	rectification::Symbol
	log_energy_source::Symbol
	log_energy_pos::Symbol
	delta_window_length::Int64
	delta_matrix::Symbol

	# spectral
	spectral_spectrum::Symbol

	# f0
	f0_method::Symbol
	f0_range::Tuple{Int64, Int64}
	median_filter_length::Int64

	# constant-q transform
	bins_octave::Int64
	freq_limits::Tuple{Float64, Float64}
	transform_type::Symbol
end

@with_kw mutable struct AudioData
	x::AbstractVector{Float64}

	# fft
	fft::AbstractArray{Float64} = []
	fft_window::Vector{Float64} = []

	# linear spectrum
	lin_spectrogram::AbstractArray{Float64} = []
	lin_frequencies::Vector{Float64} = []

	# mel_spectrum
	mel_filterbank::AbstractArray{Float64} = []
	mel_spectrogram::AbstractArray{Float64} = []
	log_mel::AbstractArray{Float64} = []

	# mfcc
	mfcc_coeffs::AbstractArray{Float64} = []
	mfcc_delta::AbstractArray{Float64} = []
	mfcc_deltadelta::AbstractArray{Float64} = []
	log_energy::Vector{Float64} = []

	# spectral
	spectral_centroid::Vector{Float64} = []
	spectral_crest::Vector{Float64} = []
	spectral_decrease::Vector{Float64} = []
	spectral_entropy::Vector{Float64} = []
	spectral_flatness::Vector{Float64} = []
	spectral_flux::Vector{Float64} = []
	spectral_kurtosis::Vector{Float64} = []
	spectral_rolloff::Vector{Float64} = []
	spectral_skewness::Vector{Float64} = []
	spectral_slope::Vector{Float64} = []
	spectral_spread::Vector{Float64} = []

	# f0
	f0::Vector{Float64} = []

	# constant-q transform
	cqt_spec::AbstractArray{Float64} = []
end

# reference:
# https://www.functionalnoise.com/pages/2023-01-31-julia-class/

mutable struct AudioObj
	setup::AudioSetup
	data::AudioData

	const get_fft::Function
	const get_lin_spec::Function
	const get_mel_spec::Function
	const get_log_mel::Function
	const get_mfcc::Function
	const get_spectrals::Function
	const get_f0::Function
	const get_cqt::Function
	const get_features::Function

	function get_fft(self::AudioObj)
		if isempty(self.data.fft)
			get_fft!(self.setup, self.data)
		end

		return self.data.fft'
	end

	function get_lin_spec(self::AudioObj)
		if isempty(self.data.lin_spectrogram)
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			lin_spectrogram!(self.setup, self.data)
		end

		return self.data.lin_spectrogram
	end

	function get_mel_spec(self::AudioObj)
		if isempty(self.data.mel_spectrogram)
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			get_mel_spec!(self.setup, self.data)
		end

		return self.data.mel_spectrogram
	end

	function get_log_mel(self::AudioObj)
		if isempty(self.data.log_mel)
			if isempty(self.data.mel_spectrogram)
				if isempty(self.data.fft)
					get_fft!(self.setup, self.data)
				end
				get_mel_spec!(self.setup, self.data)
			end
			get_log_mel!(self.setup, self.data)
		end

		return self.data.log_mel
	end

	function get_mfcc(self::AudioObj)
		if isempty(self.data.mfcc_coeffs)
			if isempty(self.data.mel_spectrogram)
				if isempty(self.data.fft)
					get_fft!(self.setup, self.data)
				end
				get_mel_spec!(self.setup, self.data)
			end
			get_mfcc!(self.setup, self.data)
			get_mfcc_deltas!(self.setup, self.data)
		end

		return hcat((self.data.mfcc_coeffs, self.data.mfcc_delta, self.data.mfcc_deltadelta)...)
	end

	function get_spectrals(self::AudioObj)
		if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			lin_spectrogram!(self.setup, self.data)
		elseif self.setup.spectral_spectrum == :mel && isempty(self.data.mel_spectrogram)
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			mel_spectrogram!(self.setup, self.data)
		end
		if isempty(self.data.spectral_centroid)
			get_spectrals!(self.setup, self.data)
		end

		return hcat(
			(
				self.data.spectral_centroid,
				self.data.spectral_crest,
				self.data.spectral_decrease,
				self.data.spectral_entropy,
				self.data.spectral_flatness,
				self.data.spectral_flux,
				self.data.spectral_kurtosis,
				self.data.spectral_rolloff,
				self.data.spectral_skewness,
				self.data.spectral_slope,
				self.data.spectral_spread,
			)...,
		)
	end

	function get_f0(self::AudioObj)
		if isempty(self.data.f0)
			get_f0!(self.setup, self.data)
		end

		return self.data.f0
	end

	function get_cqt(self::AudioObj)
		if isempty(self.data.cqt_spec)
			get_cqt_spec!(self.setup, self.data)
		end

		return self.data.cqt_spec
	end

	# --------------------------------------------------------------------------- #
	#                            get features profiles                            #
	# --------------------------------------------------------------------------- #
	function get_features(self::AudioObj, profile::Symbol = :full)
		if profile == :full
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			if isempty(self.data.mel_spectrogram)
				get_mel_spec!(self.setup, self.data)
			end
			if isempty(self.data.log_mel)
				get_log_mel!(self.setup, self.data)
			end
			if isempty(self.data.mfcc_coeffs)
				get_mfcc!(self.setup, self.data)
				get_mfcc_deltas!(self.setup, self.data)
			end
			if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
				lin_spectrogram!(self.setup, self.data)
			end
			if isempty(self.data.spectral_centroid)
				get_spectrals!(self.setup, self.data)
			end
			if isempty(self.data.f0)
				get_f0!(self.setup, self.data)
			end

			return hcat(
				(
					self.data.mel_spectrogram,
					self.data.log_mel,
					self.data.mfcc_coeffs,
					self.data.mfcc_delta,
					self.data.mfcc_deltadelta,
					self.data.spectral_centroid,
					self.data.spectral_crest,
					self.data.spectral_decrease,
					self.data.spectral_entropy,
					self.data.spectral_flatness,
					self.data.spectral_flux,
					self.data.spectral_kurtosis,
					self.data.spectral_rolloff,
					self.data.spectral_skewness,
					self.data.spectral_slope,
					self.data.spectral_spread,
					self.data.f0,
				)...,
			)

			# --------------------------------------------------------------------------- #
			#                           emotion work in progress                          #
			# --------------------------------------------------------------------------- #
		elseif profile == :emotion
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			if isempty(self.data.mel_spectrogram)
				get_mel_spec!(self.setup, self.data)
			end
			if isempty(self.data.log_mel)
				get_log_mel!(self.setup, self.data)
			end
			if isempty(self.data.mfcc_coeffs)
				get_mfcc!(self.setup, self.data)
				get_mfcc_deltas!(self.setup, self.data)
			end
			if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
				lin_spectrogram!(self.setup, self.data)
			end
			if isempty(self.data.spectral_centroid)
				get_spectrals!(self.setup, self.data)
			end
			if isempty(self.data.f0)
				get_f0!(self.setup, self.data)
			end

			return hcat(
				(
					# self.data.mel_spectrogram,
					self.data.log_mel,
					# self.data.mfcc_coeffs,
					self.data.mfcc_delta,
					self.data.mfcc_deltadelta,
					self.data.spectral_centroid,
					self.data.spectral_crest,
					self.data.spectral_decrease,
					self.data.spectral_entropy,
					self.data.spectral_flatness,
					self.data.spectral_flux,
					self.data.spectral_kurtosis,
					self.data.spectral_rolloff,
					self.data.spectral_skewness,
					self.data.spectral_slope,
					self.data.spectral_spread,
					self.data.f0,
				)...,
			)

			# elseif profile == :emotion_set_1
			# 	if isempty(self.data.fft)
			# 		get_fft!(self.setup, self.data)
			# 	end
			# 	if isempty(self.data.mel_spectrogram)
			# 		get_mel_spec!(self.setup, self.data)
			# 	end
			# 	if isempty(self.data.mfcc_coeffs)
			# 		get_mfcc!(self.setup, self.data)
			# 		# get_mfcc_deltas!(self.setup, self.data)
			# 	end
			# 	if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
			# 		lin_spectrogram!(self.setup, self.data)
			# 	end
			# 	if isempty(self.data.spectral_centroid)
			# 		get_spectrals!(self.setup, self.data)
			# 	end
			# 	# if isempty(self.data.f0)
			# 	#     get_f0!(self.setup, self.data)
			# 	# end

			# 	return vcat(
			# 		(
			# 			self.data.mel_spectrogram',
			# 			self.data.mfcc_coeffs',
			# 			# self.data.mfcc_delta',
			# 			# self.data.mfcc_deltadelta',
			# 			self.data.spectral_centroid',
			# 			self.data.spectral_crest',
			# 			self.data.spectral_decrease',
			# 			self.data.spectral_entropy',
			# 			self.data.spectral_flatness',
			# 			self.data.spectral_flux',
			# 			self.data.spectral_kurtosis',
			# 			self.data.spectral_rolloff',
			# 			self.data.spectral_skewness',
			# 			self.data.spectral_slope',
			# 			self.data.spectral_spread',
			# 			# self.data.f0'
			# 		)...,
			# 	)

			# elseif profile == :emotion_set_2
			# 	if isempty(self.data.fft)
			# 		get_fft!(self.setup, self.data)
			# 	end
			# 	if isempty(self.data.mfcc_coeffs)
			# 		get_mel_spec!(self.setup, self.data)
			# 		get_mfcc!(self.setup, self.data)
			# 		get_mfcc_deltas!(self.setup, self.data)
			# 	end
			# 	if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
			# 		lin_spectrogram!(self.setup, self.data)
			# 	end
			# 	if isempty(self.data.spectral_centroid)
			# 		get_spectrals!(self.setup, self.data)
			# 	end
			# 	if isempty(self.data.f0)
			# 		get_f0!(self.setup, self.data)
			# 	end

			# 	return vcat(
			# 		(
			# 			# self.data.mel_spectrogram',
			# 			self.data.mfcc_coeffs',
			# 			self.data.mfcc_delta',
			# 			self.data.mfcc_deltadelta',
			# 			self.data.spectral_centroid',
			# 			self.data.spectral_crest',
			# 			self.data.spectral_decrease',
			# 			self.data.spectral_entropy',
			# 			self.data.spectral_flatness',
			# 			self.data.spectral_flux',
			# 			self.data.spectral_kurtosis',
			# 			self.data.spectral_rolloff',
			# 			self.data.spectral_skewness',
			# 			self.data.spectral_slope',
			# 			self.data.spectral_spread',
			# 			self.data.f0',
			# 		)...,
			# 	)

			# elseif profile == :emotion_set_3
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			if isempty(self.data.mel_spectrogram)
				get_mel_spec!(self.setup, self.data)
			end
			if isempty(self.data.mfcc_coeffs)
				get_mfcc!(self.setup, self.data)
				get_mfcc_deltas!(self.setup, self.data)
			end
			if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
				lin_spectrogram!(self.setup, self.data)
			end
			if isempty(self.data.spectral_centroid)
				get_spectrals!(self.setup, self.data)
			end
			if isempty(self.data.f0)
				get_f0!(self.setup, self.data)
			end

			return vcat(
				(
					self.data.mel_spectrogram',
					self.data.mfcc_coeffs',
					self.data.mfcc_delta',
					self.data.mfcc_deltadelta',
					self.data.spectral_centroid',
					self.data.spectral_crest',
					self.data.spectral_decrease',
					self.data.spectral_entropy',
					self.data.spectral_flatness',
					self.data.spectral_flux',
					self.data.spectral_kurtosis',
					self.data.spectral_rolloff',
					self.data.spectral_skewness',
					self.data.spectral_slope',
					self.data.spectral_spread',
					self.data.f0',
				)...,
			)

			# --------------------------------------------------------------------------- #
			#                                 Gio Paglia                                  #
			# --------------------------------------------------------------------------- #

		elseif profile == :kdd
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			if isempty(self.data.mel_spectrogram)
				get_mel_spec!(self.setup, self.data)
			end
			if isempty(self.data.mfcc_coeffs)
				get_mfcc!(self.setup, self.data)
				# get_mfcc_deltas!(self.setup, self.data)
			end
			# if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
			# 	lin_spectrogram!(self.setup, self.data)
			# end
			# if isempty(self.data.spectral_centroid)
			# 	get_spectrals!(self.setup, self.data)
			# end
			if isempty(self.data.f0)
				self.setup.f0_range = (200, 700)
				get_f0!(self.setup, self.data)
			end

			return hcat(
				(
					# self.data.mel_spectrogram,
					self.data.mfcc_coeffs,
					# self.data.mfcc_delta,
					# self.data.mfcc_deltadelta,
					# self.data.spectral_centroid,
					# self.data.spectral_crest,
					# self.data.spectral_decrease,
					# self.data.spectral_entropy,
					# self.data.spectral_flatness,
					# self.data.spectral_flux,
					# self.data.spectral_kurtosis,
					# self.data.spectral_rolloff,
					# self.data.spectral_skewness,
					# self.data.spectral_slope,
					# self.data.spectral_spread,
					self.data.f0,
				)...,
			)
		elseif profile == :kdd2
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			if isempty(self.data.mel_spectrogram)
				get_mel_spec!(self.setup, self.data)
			end
			if isempty(self.data.log_mel)
				get_log_mel!(self.setup, self.data)
			end
			if isempty(self.data.f0)
				self.setup.f0_range = (200, 700)
				get_f0!(self.setup, self.data)
			end

			return hcat(
				(
					self.data.log_mel,
					self.data.f0,
				)...,
			)
		else
			@error("Unknown $profile profile.")
		end
	end

	function AudioObj(setup::AudioSetup, data::AudioData)
		obj = new(
			setup, data,
			() -> get_fft(obj),
			() -> get_lin_spec(obj),
			() -> get_mel_spec(obj),
			() -> get_log_mel(obj),
			() -> get_mfcc(obj),
			() -> get_spectrals(obj),
			() -> get_f0(obj),
			() -> get_cqt(obj),
			(x) -> get_features(obj, x),
		)
		#   return obj
	end
end