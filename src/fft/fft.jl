# function get_onesided_fft_range(stft_length::Int64)
# 	if mod(stft_length, 2) == 0
# 		return collect(1:Int(stft_length / 2 + 1))   # EVEN
# 	else
# 		return collect(1:Int((stft_length + 1) / 2))  # ODD
# 	end
# end # get_onesided_fft_range

# #------------------------------------------------------------------------------#
# #              fft version 1 as used in audio features extraction              #
# #------------------------------------------------------------------------------#
# function _get_stft(x::AbstractArray{Float64}, setup::AudioSetup)
# 	hop_length = setup.stft.win_length - setup.stft.overlap_length
# 	if isempty(setup.stft.win)
# 		setup.stft.win, _ = gencoswin(setup.stft.win_type[1], setup.stft.win_length, setup.stft.win_type[2])
# 	end

# 	# split in windows
# 	y = buffer(x, setup.stft.win_length, setup.stft.win_length - setup.stft.overlap_length)

# 	# apply window and take fft
# 	Z = fft(y .* setup.stft.win, (1,))

# 	# take one side
# 	logical_ossb = falses(setup.stft.stft_length)
# 	logical_ossb[get_onesided_fft_range(setup.stft.stft_length)] .= true
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

# get_stft!(setup::AudioSetup, data::AudioData) = data.stft.stft, data.log_energy = _get_stft(data.x, setup)