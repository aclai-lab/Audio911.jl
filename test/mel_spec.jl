using Test
using Audio911

using MAT

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = Audio911.load(wav_file, format=Float64)

frames = AudioFrames(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true)
stft = Stft(frames; spectrum_type=power)
mel_spec = MelSpec(
    stft;
    freq_range=(100,1000),
    win_norm=true,
    melstyle=:htk,
    nbands=26,
    fb_norm=:bandwidth,
    fb_desing=:linear,
    applylog=false
)

# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct MelSpecInfo <: Audio911.AbstractInfo
    sr            :: Int64
    freq_range    :: FreqRange
    spectrum_type :: Base.Callable
	win_norm      :: Bool
    melstyle      :: Symbol
    nbands        :: Int64
    fb_norm       :: Symbol
    fb_desing     :: Symbol
    applylog      :: Bool
end

# ---------------------------------------------------------------------------- #
#                           mel spectrogram struct                             #
# ---------------------------------------------------------------------------- #
struct MelSpec{F,T} <: Audio911.AbstractSpectrogram
    spec :: Matrix{T}
    freq :: Audio911.StepRangeLen
    info :: MelSpecInfo
    # fbank :: FbInfo

    function MelSpec{F}(
        spec :: Matrix{T},
        freq :: StepRangeLen,
        info :: MelSpecInfo
    ) where {F,T}
        new{F,T}(spec, freq, info)
    end
end

function LinSpec(
	stft       :: Stft;
	freq_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr>>1),
	win_norm   :: Bool=false
)::LinSpec
	spec = get_spec(stft)
	freq = get_freq(stft)

	sr            = get_sr(stft)
	stft_size     = get_stft_size(stft)
	spectrum_type = get_spec_type(stft)
	window        = get_window(stft)

	freq_range isa Tuple && (freq_range = FreqRange(first(freq_range), last(freq_range)))

	if freq_range != FreqRange(0, sr >> 1)
		bin_low, bin_high = get_freq_range(freq_range, stft_size, sr)
		spec = @views spec[bin_low:bin_high, :]
		freq = freq[bin_low:bin_high]
	end

	win_norm_func = eval(Symbol("win" * string(spectrum_type)))
	win_norm && (spec = win_norm_func(spec, window))

	info = LinSpecInfo(sr, freq_range, spectrum_type, win_norm)

	return LinSpec{typeof(stft)}(spec .* 2, freq, info)
end

# ---------------------------------------------------------------------------- #
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
function MelSpec(
	stft::Stft;
	freq_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr>>1),
	win_norm   :: Bool=false,
	melstyle   :: Symbol=:htk, 	       # :htk, :slaney
    nbands     :: Int64=26,
	fb_norm    :: Symbol=:bandwidth,   # :bandwidth, :area, :none
	fb_desing  :: Symbol=:linear,
	# frequency_scale  :: Symbol=:mel, # TODO :mel, :bark, :erb
    applylog   :: Bool=false
# )::MelSpec
)
    data = get_stft(stft)
    win_size = get_info(stft).win_size
    win_step = get_info(stft).win_step

	fb, mel_freq = design_fb(
        stft;
        mel_bands,
        mel_style,
        fb_design_domain,
        fb_norm,
        frequency_scale,
        st_peak_range
    )

	num_hops = size(data, 2)

	# apply fb
	# if (setup.spectrum_type == :power)
	mel_spec = reshape(fb * data, mel_bands, num_hops)
	# else
	#     #TODO
	#     error("magnitude not yet implemented.")
	# end

	# mel_spec = transpose(mel_spec)

    return mel_spec, fb, mel_freq
end