# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct LinSpecSetup <: AbstractSetup
    sr            :: Int64
    freqrange     :: FreqRange
    spectrum_type :: Base.Callable
	win_norm      :: Bool
end

# ---------------------------------------------------------------------------- #
#                          linear spectrogram struct                           #
# ---------------------------------------------------------------------------- #
struct LinSpec{F,T} <: AbstractSpectrogram
    spec :: Matrix{T}
    freq :: StepRangeLen
    info :: LinSpecSetup

    function LinSpec{F}(
        spec :: Matrix{T},
        freq :: StepRangeLen,
        info :: LinSpecSetup
    ) where {F,T}
        new{F,T}(spec, freq, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                  utilities                                   #
# ---------------------------------------------------------------------------- #
function get_freq_range(
    freqrange :: FreqRange,
    stft_size  :: Int64,
    sr         :: Int64
)
    # convert frequencies to bin indices
    bin_low  = cld(get_low(freqrange) * stft_size, sr) + 1
    bin_high = fld(get_hi(freqrange)  * stft_size, sr) + 1

    return bin_low, bin_high
end

# ---------------------------------------------------------------------------- #
#                                 get lin spec                                 #
# ---------------------------------------------------------------------------- #
function LinSpec(
	stft       :: Stft;
	freqrange :: FreqRange=(0, get_setup(frames).sr>>1),
	win_norm   :: Bool=false
)::LinSpec
	spec = get_data(stft)
	freq = get_freq(stft)

	sr            = get_sr(stft)
	stft_size     = get_nfft(stft)
	spectrum_type = get_spectrum(stft)
	window        = get_window(stft)

	if freqrange != (0, sr >> 1)
		bin_low, bin_high = get_freq_range(freqrange, stft_size, sr)
		spec = @views spec[bin_low:bin_high, :]
		freq = freq[bin_low:bin_high]
	end

	win_norm_func = eval(Symbol("win" * string(spectrum_type)))
	win_norm && (spec = win_norm_func(spec, window))

	info = LinSpecSetup(sr, freqrange, spectrum_type, win_norm)

	return LinSpec{typeof(stft)}(spec .* 2, freq, info)
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
Base.eltype(::LinSpec{T}) where T = T

get_data(s::LinSpec)  = s.spec'
get_freq(s::LinSpec)  = s.freq
get_setup(s::LinSpec) = s.info
