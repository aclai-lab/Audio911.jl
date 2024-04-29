"""
	Audio911 signal data structures

	uses package Parameters for @with_kw mutable struct

	AudioSetup stores all datas that has to be shared in Audio911 module
	AudioData stores all results from signal analysis
"""
@with_kw mutable struct AudioSetup
	sr::Int64

	# fft
	fft_length::Int64 = 0
	window_type::Tuple{Symbol, Symbol} = (:hann, :periodic)
	window_length::Int64 = 0
	overlap_length::Int64 = 0
	window_norm::Bool = true

	# spectrum
	frequency_range::Tuple{Int64, Int64} = (0, 0)
	lin_frequencies::Vector{Float64} = []
	band_edges::AbstractVector{AbstractFloat} = []
	spectrum_type::Symbol = :power # :power, :magnitude

	# mel
	mel_style::Symbol = :htk # :htk, :slaney
	mel_bands::Int64 = 26
	mel_frequencies::Vector{Float64} = []
	filterbank_design_domain::Symbol = :linear
	filterbank_normalization::Symbol = :bandwidth # :bandwidth, :area, :none
	frequency_scale::Symbol = :mel # TODO :mel, :bark, :erb

	# mfcc
	num_coeffs::Int64 = 13
	normalization_type::Symbol = :dithered # :standard, :dithered
	rectification::Symbol = :log # :log, :cubic_root
	log_energy_source::Symbol = :standard # :standard (after windowing), :mfcc
	log_energy_pos::Symbol = :none #:append, :replace, :none
	delta_window_length::Int64 = 9
	delta_matrix::Symbol = :transposed # :standard, :transposed

	# spectral
	spectral_spectrum::Symbol = :lin # :lin, :mel

	# f0
	f0_method::Symbol = :nfc
	f0_range::Tuple{Int64, Int64} = (50, 400)
	median_filter_length::Int64 = 1
end

@with_kw mutable struct AudioData
	x::AbstractVector{Float64} = []

	# fft
	fft::AbstractArray{Float64} = []
	fft_window::Vector{Float64} = []

	# linear spectrum
	lin_spectrogram::AbstractArray{Float64} = []
	lin_frequencies::Vector{Float64} = []

	# mel_spectrum
	mel_filterbank::AbstractArray{Float64} = []
	mel_spectrogram::AbstractArray{Float64} = []

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
end

# reference:
# https://www.functionalnoise.com/pages/2023-01-31-julia-class/

# TODO metti a posto gli output dei metodi

mutable struct AudioObj
	setup::AudioSetup
	data::AudioData

	const get_fft::Function
	const get_lin_spec::Function
	const get_mel_spec::Function
	const get_mfcc::Function
	const get_spectrals::Function
	const get_f0::Function
	const get_features::Function

	function get_fft(self::AudioObj)
		if isempty(self.data.fft)
			get_fft!(self.setup, self.data)
		end

		return self.data.fft
	end

	function get_lin_spec(self::AudioObj)
		if isempty(self.data.lin_spectrogram)
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			lin_spectrogram!(self.setup, self.data)
		end

		return self.data.lin_spectrogram, self.setup.lin_frequencies
	end

	function get_mel_spec(self::AudioObj)
		if isempty(self.data.mel_spectrogram)
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			get_mel_spec!(self.setup, self.data)
		end

		return self.data.mel_spectrogram, self.setup.mel_frequencies
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

		return self.data.mfcc_coeffs, self.data.mfcc_delta, self.data.mfcc_deltadelta
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

		return [
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
		]
	end

	function get_f0(self::AudioObj)
		if isempty(self.data.f0)
			get_f0!(self.setup, self.data)
		end

		return self.data.fft
	end

	# --------------------------------------------------------------------------- #
	#                            get features profiles                            #
	# --------------------------------------------------------------------------- #
	function get_features(self::AudioObj; profile::Symbol)
		if profile == :full
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
			#                           emotion work in progress                          #
			# --------------------------------------------------------------------------- #
		elseif profile == :emotion_set_1
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
			if self.setup.spectral_spectrum == :lin && isempty(self.data.lin_spectrogram)
				lin_spectrogram!(self.setup, self.data)
			end
			if isempty(self.data.spectral_centroid)
				get_spectrals!(self.setup, self.data)
			end
			# if isempty(self.data.f0)
			#     get_f0!(self.setup, self.data)
			# end

			return vcat(
				(
					self.data.mel_spectrogram',
					self.data.mfcc_coeffs',
					# self.data.mfcc_delta',
					# self.data.mfcc_deltadelta',
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
					# self.data.f0'
				)...,
			)

		elseif profile == :emotion_set_2
			if isempty(self.data.fft)
				get_fft!(self.setup, self.data)
			end
			if isempty(self.data.mfcc_coeffs)
				get_mel_spec!(self.setup, self.data)
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
					# self.data.mel_spectrogram',
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

		elseif profile == :emotion_set_3
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
					self.data.f0'
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
				get_f0!(self.setup, self.data)
			end

			return vcat(
				(
					# self.data.mel_spectrogram',
					self.data.mfcc_coeffs',
					# self.data.mfcc_delta',
					# self.data.mfcc_deltadelta',
					# self.data.spectral_centroid',
					# self.data.spectral_crest',
					# self.data.spectral_decrease',
					# self.data.spectral_entropy',
					# self.data.spectral_flatness',
					# self.data.spectral_flux',
					# self.data.spectral_kurtosis',
					# self.data.spectral_rolloff',
					# self.data.spectral_skewness',
					# self.data.spectral_slope',
					# self.data.spectral_spread',
					self.data.f0',
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
			() -> get_mfcc(obj),
			() -> get_spectrals(obj),
			() -> get_f0(obj),
			(x) -> get_features(obj; profile = x),
		)
		#   return obj
	end
end

mutable struct MyClass
	myInt::Int

	# we have these `const` fields since Julia 1.8
	const print_int::Function
	const set_int!::Function

	function print_int(self::MyClass)
		println("hello, I have myInt: $(self.myInt)")
	end

	function set_int!(self::MyClass, new_int::Int)
		self.myInt = new_int
		return self
	end

	function MyClass(int::Int)
		obj = new(
			int,
			() -> print_int(obj),
			(new_int,) -> set_int!(obj, new_int),
		)
		return obj
	end
end

################################################################################
#                            new data structures                               #
#                                    TODO                                      #
################################################################################

# mutable struct SignalSetup
#     sr::Int

#     # fft
#     fft_length::Int
#     window_type::Tuple{Symbol, Symbol} # (:hann, :periodic)
#     window_length::Int
#     overlap_length::Int
#     window_norm::Bool

#     # spectrum
#     frequency_range::Tuple{Int64, Int64}
#     lin_frequencies::Vector{AbstractFloat}
#     band_edges::AbstractVector{AbstractFloat}
#     spectrum_type::Symbol

#     # # mel
#     # mel_style::Symbol = :htk # :htk, :slaney
#     # mel_bands::Int64 = 26
#     # mel_frequencies::Vector{Float64} = []
#     # filterbank_design_domain::Symbol = :linear
#     # filterbank_normalization::Symbol = :bandwidth # :bandwidth, :area, :none
#     # frequency_scale::Symbol = :mel # TODO :mel, :bark, :erb

#     # # mfcc
#     # num_coeffs::Int64 = 13
#     # normalization_type::Symbol = :standard # :standard, :dithered
#     # rectification::Symbol = :log # :log, :cubic_root
#     # log_energy_source::Symbol = :standard # :standard (after windowing), :mfcc
#     # log_energy_pos::Symbol = :append #:append, :replace, :none
#     # delta_window_length::Int64 = 9
#     # delta_matrix::Symbol = :standard # :standard, :transposed

#     # # spectral
#     # spectral_spectrum::Symbol = :linear # :linear, :mel

#     function SignalSetup(;
#         sr::Int,
#         fft_length::Int,
#         window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
#         window_length::Int = 0,
#         overlap_length::Int = 0,
#         window_norm::Bool = false,
#         # spectrum
#         frequency_range::Tuple{Int64, Int64} = (0, 0),
#         lin_frequencies::Vector{AbstractFloat} = [],
#         band_edges::Vector{AbstractFloat} = [],
#         spectrum_type::Symbol=:power,
#     )
#     if window_type[1] ∉ (:hann, :hamming, :blackman, :flattopwin, :rect)
#         error("Unknown window_type $window_type[1].")
#     end
#     if window_type[2] ∉ (:periodic, :symmetric)
#         error("window_type second parameter must be :periodic or :symmetric.")
#     end

#     if window_length == 0
#         window_length = fft_length
#     elseif window_length < fft_length
#         error("window_length can't be smaller than fft_length.")
#     end

#     if overlap_length == 0
#         overlap_length=Int(round(FFTLength * 0.500))
#     elseif overlap_length > window_length
#         error("overlap_length can't be greater than window_length.")
#     end

#     # if isempty(frequency_range)

#     if spectrum_type ∉ (:power, :magnitude)
#         error("spectrum_type parameter must be symbol, :power or :magnitude.")
#     end

#         new(
#             sr,

#             # fft
#             fft_length,
#             window_type,
#             window_length,
#             overlap_length,
#             window_norm,

#             # spectrum
#             frequency_range,
#             lin_frequencies,
#             band_edges,
#             spectrum_type,

#             # # mel
#             # mel_style::Symbol = :htk # :htk, :slaney
#             # mel_bands::Int64 = 26
#             # mel_frequencies::Vector{Float64} = []
#             # filterbank_design_domain::Symbol = :linear
#             # filterbank_normalization::Symbol = :bandwidth # :bandwidth, :area, :none
#             # frequency_scale::Symbol = :mel # TODO :mel, :bark, :erb

#             # # mfcc
#             # num_coeffs::Int64 = 13
#             # normalization_type::Symbol = :standard # :standard, :dithered
#             # rectification::Symbol = :log # :log, :cubic_root
#             # log_energy_source::Symbol = :standard # :standard (after windowing), :mfcc
#             # log_energy_pos::Symbol = :append #:append, :replace, :none
#             # delta_window_length::Int64 = 9
#             # delta_matrix::Symbol = :standard # :standard, :transposed

#             # # spectral
#             # spectral_spectrum::Symbol = :linear # :linear, :mel
#         )
#     end
# end