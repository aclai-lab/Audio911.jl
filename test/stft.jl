#------------------------------------------------------------------------------#
#                                  utilities                                   #
#------------------------------------------------------------------------------#
function _trim_freq_range(
    x::Union{AbstractArray{Float64}, NamedArray{Float64}};
    sr::Int64,
    win_length::Int64,
    freq_range::Tuple{Int64, Int64}
)
    # trim to desired range
    bin_low = ceil(Int, freq_range[1] * win_length / sr + 1)
    bin_high = floor(Int, freq_range[2] * win_length / sr + 1)
    bins = bin_low:bin_high

    # create frequency vector
    freq_vec = floor.(Int, (sr / win_length) .* (collect(bins) .- 1))
    # shift final bin if fftLength is odd and the final range is full to sr/2.
    if (rem(win_length, 2) == 1 && bin_high == floor(win_length / 2 + 1))
        freq_vec[end] = floor(Int, sr * (win_length - 1) / (2 * win_length))
    end

    return x[bins, :], freq_vec, bins
end

function _trim_freq_range(x::Union{AbstractArray{Float64}, NamedArray{Float64}}, s::StftSetup)
    _trim_freq_range(x; sr=s.sr, win_length=s.stft_length, freq_range=s.freq_range)
end

function _trim_freq_range(x::Union{AbstractArray{Float64}, NamedArray{Float64}}, a::AudioSetup)
    _trim_freq_range(x; sr=a.sr, win_length=a.stft_length, freq_range=a.freq_range)  
end

#------------------------------------------------------------------------------#
#                                     stft                                     #
#------------------------------------------------------------------------------#
"""
Note on choose win_length and fft_length values

When you apply a Fast Fourier Transform (FFT), you're transforming your data from the time domain to the frequency domain.
The length of the FFT, here called fft_length, determines the resolution of the frequency domain output.
If the window length (the length of the data segment you're analyzing) is greater than the FFT length, 
you have to somehow reduce the amount of data to fit into the FFT. One common way to do this is by "wrapping" the data.

Wrapping involves taking the part of the data that "overflows" the FFT length and adding it back to the beginning 
of the data segment. This is equivalent to assuming that the data is periodic and continues to repeat.
However, this can introduce discontinuities at the wrap-around point, 
which can lead to spectral leakage (distortion in the frequency domain representation). 
To mitigate this, it's common to apply a window function to the data before performing the FFT.

It's important to choose an appropriate FFT length for your specific application. 
If you're dealing with non-periodic signals or you want to avoid the potential issues associated with wrapping, 
it might be better to choose an FFT length that's equal to or larger than your window length.

If the FFT window is larger than the window, the audio data will be zero-padded to match the size of the FFT window.
This zero-padding in the time domain results in an interpolation in the frequency domain, 
which can provide a more detailed view of the spectral content of the signal.
"""
function _get_stft(
        wframes::Union{AbstractArray{Float64}, NamedArray{Float64}};
        sr::Int64,
        win::AbstractVector{Float64},
        win_length::Int64,
        stft_length::Int64,
        col_offsets::AbstractVector{Int64},
        spec_norm::Symbol,
        apply_log::Bool
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
    stft_freq = floor.(Int, (sr / stft_length) * (0:size(stft_spec, 1) .- 1))

    spectrum_funcs = Dict(
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :winpower => x -> real.((x .* conj.(x))) / sum(win)^2,
        :winmagnitude => x -> abs.(x) / sum(win),
	)
    # check if spectrum_type is valid
    @assert haskey(spectrum_funcs, spec_norm) "Unknown spectrum_type: $spec_norm."
    stft_spec = spectrum_funcs[spec_norm](stft_spec)

    if apply_log
        stft_spec = log10.(stft_spec .+ 1e-10)
    end

    NamedArray(stft_spec; names=(collect(stft_freq)u"Hz", collect(col_offsets)))
end

function _get_stft(wframes::AbstractArray{Float64}, s::StftSetup)
    _get_stft(wframes; sr=s.sr, win=s.win, win_length=s.win_length, stft_length=s.stft_length, col_offsets=s.col_offsets, spec_norm=s.spec_norm, apply_log=s.apply_log)
end

function _get_stft(wframes::AbstractArray{Float64}, a::AudioSetup)
    _get_stft(wframes; sr=a.sr, win=a.win, win_length=a.win_length, stft_length=a.stft_length, col_offsets=a.col_offsets, spec_norm=a.spec_norm, apply_log=a.apply_log)
end

function get_stft(
        x::AbstractVector{Float64},
        sr::Int64;
        stft_length::Int64 = sr <= 8000 ? 256 : 512,
        win_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
        win_length::Int64 = stft_length,
        overlap_length::Int64 = round(Int, win_length / 2),
        freq_range::Tuple{Int64, Int64} = (0, floor(Int, sr / 2)),
        spec_norm::Symbol = :power, # :none, :power, :magnitude, :winpower, :winmagnitude
        apply_log::Bool = false
)
    # create stft setup object
    stft_setup = StftSetup(
        sr=sr,
        stft_length=stft_length,
        win_type=win_type,
        win_length=win_length,
        overlap_length=overlap_length,
        freq_range=freq_range,
        spec_norm=spec_norm,
        apply_log=apply_log
    )

    # get partitioned signal into frames, and window coefficients
    _, stft_setup.win, wframes, _, stft_setup.col_offsets = _get_frames(x, stft_setup)

    # if stft_length > win_length
    #     @info("FFT window size is greater than actual window size. Zero-padding will be applied.")
    # end

    # get stft
    stft_spec = _get_stft(wframes, stft_setup)

    stft_spec, _, _ = _trim_freq_range(stft_spec, stft_setup)

    # adjust bins
    # if freq_range[1] == 0
    #     stft_spec[1, :] = stft_spec[1, :] .* 0.5
    # end
    # if freq_range[2] == floor(Int, sr / 2)
    #     stft_spec[end, :] = stft_spec[end, :] .* 0.5
    # end

    return stft_spec
end

function get_stft(x::AbstractVector{<:AbstractFloat}, sr::Int64; kwargs...)
    get_stft(Float64.(x), sr; kwargs...)
end

function get_stft!(a::AudioObj)
	_get_stft(a.data.frames .* a.setup.window, a.setup)
end


# -------------------------------------------------------------------------- #
#                                   debug                                    #
# -------------------------------------------------------------------------- #
using Revise, Plots
using Audio911, Parameters, FFTW, DSP, StatsBase, NaNStatistics
using Unitful, NamedArrays

TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
# TESTFILE = "104_1b1_Al_sc_Litt3200_4.wav"
wavfile = joinpath(TESTPATH, TESTFILE)

sr_src = 16000
x, sr = load_audio(wavfile, sr_src)
x = Float64.(x)

include("/home/paso/.julia/dev/Audio911.jl/src/windowing/windows.jl")
include("/home/paso/.julia/dev/Audio911.jl/src/windowing/windowing.jl")
# extra parameters:
stft_length = 1024
win_type = (:hamming, :symmetric)
win_length = 512
overlap_length = 256
spec_norm = :power
apply_log = false
freq_range = (0, 8000)

stft_spec, stft_freq = get_stft(x, sr; stft_length=stft_length, win_type=win_type, win_length=win_length, overlap_length=overlap_length, freq_range=freq_range, spec_norm=spec_norm, apply_log=apply_log)