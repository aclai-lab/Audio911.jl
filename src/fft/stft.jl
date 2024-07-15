# Note on choose win_length and stft_length values

# When you apply a Fast Fourier Transform (FFT), you're transforming your data from the time domain to the frequency domain.
# The length of the FFT, here called stft_length, determines the resolution of the frequency domain output.
# If the window length (the length of the data segment you're analyzing) is greater than the FFT length, 
# you have to somehow reduce the amount of data to fit into the FFT. One common way to do this is by "wrapping" the data.

# Wrapping involves taking the part of the data that "overflows" the FFT length and adding it back to the beginning 
# of the data segment. This is equivalent to assuming that the data is periodic and continues to repeat.
# However, this can introduce discontinuities at the wrap-around point, 
# which can lead to spectral leakage (distortion in the frequency domain representation). 
# To mitigate this, it's common to apply a window function to the data before performing the FFT.

# It's important to choose an appropriate FFT length for your specific application. 
# If you're dealing with non-periodic signals or you want to avoid the potential issues associated with wrapping, 
# it might be better to choose an FFT length that's equal to or larger than your window length.

# If the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
# This zero-padding in the time domain results in an interpolation in the frequency domain, 
# which can provide a more detailed view of the spectral content of the signal.

# ---------------------------------------------------------------------------- #
#                                     stft                                     #
# ---------------------------------------------------------------------------- #
function _get_stft(
        wframes::AbstractArray{Float64};
        sr::Int64,
        win_length::Int64,
        stft_length::Int64,
        spec_norm::Symbol
)
    if win_length < stft_length
        wframes = vcat(wframes, zeros(Float64, stft_length - win_length, size(wframes, 2)))
    elseif win_length > stft_length
        @error("FFT window size smaller than actual window size is highly discuraged. Consider to review your setup.")
    end
    stft_spec = fft(wframes, (1,))

    # take one side
    if mod(stft_length, 2) == 0
        one_side = [1:Int(stft_length / 2 + 1), :]   # even
    else
        one_side = [1:Int((stft_length + 1) / 2), :]  # odd
    end
    stft_spec = stft_spec[one_side...]

    # get frequency vector
    stft_freq = (sr / stft_length) * (0:(size(stft_spec, 1) .- 1))

    # normalize
    norm_funcs = Dict(
        :none => x -> x,
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => sqrt.(real.((x .* conj.(x)))),
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, spec_norm) "Unknown spectrum_type: $spec_norm."
    stft_spec = norm_funcs[spec_norm](stft_spec)

    return stft_spec, stft_freq
end

function _get_stft!(
    rack::AudioRack;
    stft_length::Int64 = rack.audio.sr <= 8000 ? 256 : 512,
    win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
    win_length::Int64 = stft_length, 				   	 # standard setting: round(Int, 0.03 * sr)
	overlap_length::Int64 = round(Int, stft_length / 2), # standard setting: round(Int, 0.02 * sr)
	norm::Symbol = :power, # :none, :power, :magnitude, :pow2mag
)
    frames, win, wframes, _, _ = _get_frames(
        rack.audio.data,
        win_type=win_type,
        win_length=win_length,
        overlap_length=overlap_length
    )

    stft, freq = _get_stft(
        wframes,
        sr = rack.audio.sr,
        win_length = win_length,
        stft_length = stft_length,
        spec_norm = norm
    ) 

    rack.stft = Stft(
        stft_length,
        win,
        win_type,
        win_length,
        overlap_length,
        norm,
        frames,
        stft,
        freq
    )
end

# ---------------------------------------------------------------------------- #
#                                  callings                                    #
# ---------------------------------------------------------------------------- #
function get_stft!(
    rack::AudioRack,
    audio::Audio;
    kwargs...
)
    _get_stft!(rack; kwargs...)
end

# ---------------------------------------------------------------------------- #
#                            linear spectrogram                                #
# ---------------------------------------------------------------------------- #
no_zero =  x -> x == 0 ? floatmin(Float64) : x

function _get_lin(
    x::AbstractArray{Float64},
    win::AbstractVector{Float64},
    x_freq::StepRangeLen{Float64};
    freq_range::Tuple{Int64, Int64},
    win_norm::Symbol,
    db_scale::Bool
)
    # trim to desired range
    x_range = findall(freq_range[1] .<= x_freq .<= freq_range[2])
    lin_spec, lin_freq = x[x_range, :], x_freq[x_range]

    # normalize
    norm_funcs = Dict(
        :none => x -> x,
        :power => x -> x / sum(win)^2,
        :magnitude => x / sum(win),
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, win_norm) "Unknown spectrum_type: $spec_norm."
    lin_spec = norm_funcs[win_norm](x * 2)

    # scale to log
    if db_scale
        lin_spec = log10.(no_zero.(lin_spec))
    end

    return lin_spec, lin_freq
end

function _get_lin!(
    rack::AudioRack;
    freq_range::Tuple{Int64, Int64} = (0, floor(Int, rack.audio.sr / 2)),
    win_norm::Symbol = :power,
    db_scale::Bool = false
)
    spec, freq = _get_lin(
        rack.stft.stft,
        rack.stft.win,
        rack.stft.freq,
        freq_range=freq_range,
        win_norm=win_norm,
        db_scale=db_scale
    )

    rack.lin = LinSpec(
        win_norm,
        db_scale,
        spec,
        freq
    )
end

# ---------------------------------------------------------------------------- #
#                                  callings                                    #
# ---------------------------------------------------------------------------- #
function get_lin!(
    rack::AudioRack,
    stft::Stft;
    kwargs...
)
    _get_lin!(rack; kwargs...)
end

# TODO: calling for wavelets