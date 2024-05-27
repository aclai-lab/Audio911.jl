function get_onesided_fft_range(fft_length::Int64)
	if mod(fft_length, 2) == 0
		return collect(1:Int(fft_length / 2 + 1))   # EVEN
	else
		return collect(1:Int((fft_length + 1) / 2))  # ODD
	end
end # get_onesided_fft_range

#------------------------------------------------------------------------------#
#              fft version 1 as used in audio features extraction              #
#------------------------------------------------------------------------------#
function _get_fft(x::AbstractArray{Float64}, setup::AudioSetup)
	hop_length = setup.win_length - setup.overlap_length
	if isempty(setup.win)
		setup.win, _ = gencoswin(setup.win_type[1], setup.win_length, setup.win_type[2])
	end

	# split in windows
	y = buffer(x, setup.win_length, setup.win_length - setup.overlap_length)

	# apply window and take fft
	Z = fft(y .* setup.win, (1,))

	# take one side
	logical_ossb = falses(setup.fft_length)
	logical_ossb[get_onesided_fft_range(setup.fft_length)] .= true
	Z = Z[logical_ossb, :]

	# log energy
	# reference: ETSI ES 201 108 V1.1.2 (2000-04)
	# https://www.3gpp.org/ftp/tsg_sa/TSG_SA/TSGS_13/Docs/PDF/SP-010566.pdf
	if setup.log_energy_pos != :none && setup.log_energy_source == :standard
		log_energy = sum(eachrow(y .^ 2))

		if setup.normalization_type == :standard
			log_energy[log_energy.==0] .= floatmin(Float64)
		elseif setup.normalization_type == :dithered
			log_energy[log_energy.<1e-8] .= 1e-8
		end
		log_energy = log.(log_energy)
	else
		log_energy = Vector{Float64}()
	end

	if setup.spectrum_type == :power
		real(Z .* conj(Z)), log_energy
	elseif setup.spectrum_type == :magnitude
		abs.(Z), log_energy
	else
		error("Unknown spectrum type: $(setup.spectrum_type)")
	end
end

get_fft!(setup::AudioSetup, data::AudioData) = data.fft, data.log_energy = _get_fft(data.x, setup)