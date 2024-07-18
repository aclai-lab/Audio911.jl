# ---------------------------------------------------------------------------- #
#                        continuous wavelets transform                         #
# ---------------------------------------------------------------------------- #

function _get_cwt(x::AbstractVector{Float64}, fbank::AbstractArray{Float64})
    # fourier transform of input
    # obtain the CWT in the Fourier domain
    # invert to obtain wavelet coefficients
    ifft(fft(x)[1:size(fbank, 2)]' .* fbank, 2)

    # fft impleementation in matlab differs from julia
    # this code should be used for strong compatibility
    # but keep in mind that reversing the sign in imm part
    # will reverse the spectrogram
    # xposdft = fft(x)
    # cfsposdft = xposdft' .* fbank
    # cfsposdft_rev_imm = map(z -> real(z) - imag(z) * im, cfsposdft)
    # reverse(ifft(cfsposdft_rev_imm, 2))
end

function _get_cwt!(
    rack::AudioRack;
    wavelet::Symbol = :morse,
    morse_params::Tuple{Int64, Int64} = (3,60),
    vpo::Int64 = 10,
    freq_range::Tuple{Int64, Int64} = (0, floor(Int, rack.audio.sr / 2))
)
    if isnothing(rack.cwt_fb)
        get_cwt_fb!(rack, rack.audio; wavelet, morse_params, vpo, freq_range)
    end

    cwt = _get_cwt(rack.audio.data, rack.cwt_fb.fbank)

    rack.cwt = Cwt(reverse(cwt))
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