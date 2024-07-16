# ---------------------------------------------------------------------------- #
#                        continuous wavelets transform                         #
# ---------------------------------------------------------------------------- #

function _get_cwt(
        x::AbstractVector{Float64},
        x_length::Int64,
        signal_pad::Int64,
        fbank::AbstractArray{Float64}
)
    if (signal_pad > 0)
        x = vcat(reverse(x[1:signal_pad]), x, x[end:-1:(end - signal_pad + 1)])
    end

    # fourier transform of input
    xposdft = fft(x)
    # obtain the CWT in the Fourier domain
    cfsposdft = xposdft' .* fbank
    # reverse imm sign to conform matlab
    cfsposdft_rev_imm = map(z -> real(z) - imag(z) * im, cfsposdft)
    # invert to obtain wavelet coefficients
    cfs = ifft(cfsposdft_rev_imm, 2)

    if (signal_pad > 0)
        return cfs[:, signal_pad + 1:signal_pad + x_length]
    else
        return cfs
    end
end

function _get_cwt!(
    rack::AudioRack;
    wavelet::Symbol = :morse,
    morse_params::Tuple{Int64, Int64} = (3,20),
    vpo::Int64 = 10,
    boundary::Symbol = :reflection
)
    if isnothing(rack.cwt_fb)
        get_cwt_fb!(rack, rack.audio; wavelet, morse_params, vpo, boundary)
    end

    cwt = _get_cwt(rack.audio.data, size(rack.audio.data, 1), rack.cwt_fb.signal_pad, rack.cwt_fb.fbank)

    rack.cwt = Cwt(reverse(cwt, dims=1))
end

# ---------------------------------------------------------------------------- #
#                                  callings                                    #
# ---------------------------------------------------------------------------- #
function get_cwt!(
    rack::AudioRack,
    audio::Audio;
    kwargs...
)
    _get_cwt!(rack; kwargs...)
end