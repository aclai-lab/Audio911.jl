# include("../signalDataStructure.jl")

function get_lin_norm_factor(spectrum_type::Symbol, fft_window::Vector{Float64})
    if spectrum_type == :power
        return 1 / (sum(fft_window)^2) * 2
    elseif spectrum_type == :magnitude
        return 1 / sum(fft_window) * 2
    else # :none
        return 2
    end
end

################################################################################
#                                    main                                      #
################################################################################
function lin_spectrogram!(
    setup::AudioSetup,
    data::AudioData,
)
    # trim to desired range
    bin_low = Int(ceil(setup.frequency_range[1] * setup.fft_length / setup.sr + 1))
    bin_high = Int(floor(setup.frequency_range[2] * setup.fft_length / setup.sr + 1))
    bins = bin_low:bin_high

    # take one side
    bins_logical = falses(length(get_onesided_fft_range(setup.fft_length)))
    bins_logical[bins] .= true

    # create frequency vector
    w = ((setup.sr / setup.fft_length) .* (collect(bins) .- 1))
    # shift final bin if fftLength is odd and the final range is full to fs/2.
    if (rem(setup.fft_length, 2) == 1 && bin_high == floor(setup.fft_length / 2 + 1))
        w[end] = setup.sr * (setup.fft_length - 1) / (2 * setup.fft_length)
    end
    data.lin_frequencies = w[:]

    # aggiusta bordi banda, preso da audio features extraction
    adjust_first_bin = bins[1] == 1
    adjust_last_bin = (bins[end] == floor(setup.fft_length / 2 + 1)) && (rem(setup.fft_length, 2) == 0)

    if adjust_first_bin || adjust_last_bin
        adjustment = ones(length(bins))
        adjust_bins = true
        if adjust_first_bin
            adjustment[1] = 0.5
        end
        if adjust_last_bin
            adjustment[end] = 0.5
        end
    else
        adjust_bins = false
    end

    if (bins[1] == 1 && bins[end] == floor(setup.fft_length / 2 + 1))
        full_spectrum = true
    else
        full_spectrum = false
    end

    # calculate linear spectrum and normalization
    linear_norm_factor = get_lin_norm_factor(setup.spectrum_type, setup.win)
    full_spectrum ? data.lin_spectrogram = data.fft * linear_norm_factor : data.lin_spectrogram = data.fft[bins_logical, :] * linear_norm_factor

    if adjust_bins
        data.lin_spectrogram = data.lin_spectrogram .* adjustment
    end

    # linear_fc = (setup.sr ./ setup.fft_length) .* (bins .- 1)

    # # if the final bin is Nyquist, and FFTLength is even, halve it.
    # if (setup.bins[end] == floor(setup.fft_length / 2 + 1) && rem(setup.fft_length, 2) != 0)
    #     linear_fc[end] = setup.sr * (setup.fft_length - 1) / (2 * setup.fft_length)
    # end

    data.lin_spectrogram = transpose(data.lin_spectrogram)  
end # lin_spectrogram