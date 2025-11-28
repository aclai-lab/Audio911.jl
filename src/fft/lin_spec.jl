# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct LinSpecSetup <: AbstractSetup
    sr            :: Int64
    freq_range    :: FreqRange
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
#                           spectrum normalizations                            #
# ---------------------------------------------------------------------------- #
winpower(f, w)     = f / sum(w).^2
winmagnitude(f, w) = f / sum(w)

# ---------------------------------------------------------------------------- #
#                                  utilities                                   #
# ---------------------------------------------------------------------------- #
function get_freq_range(
    freq_range :: FreqRange,
    stft_size  :: Int64,
    sr         :: Int64
)
    # convert frequencies to bin indices
    bin_low  = cld(get_low(freq_range) * stft_size, sr) + 1
    bin_high = fld(get_hi(freq_range)  * stft_size, sr) + 1

    return bin_low, bin_high
end

# ---------------------------------------------------------------------------- #
#                                 get lin spec                                 #
# ---------------------------------------------------------------------------- #
function LinSpec(
	stft       :: Stft;
	freq_range :: FreqRange=(0, get_info(frames).sr>>1),
	win_norm   :: Bool=false
)::LinSpec
	spec = get_data(stft)
	freq = get_freq(stft)

	sr            = get_sr(stft)
	stft_size     = get_nfft(stft)
	spectrum_type = get_spectype(stft)
	window        = get_window(stft)

	if freq_range != (0, sr >> 1)
		bin_low, bin_high = get_freq_range(freq_range, stft_size, sr)
		spec = @views spec[bin_low:bin_high, :]
		freq = freq[bin_low:bin_high]
	end

	win_norm_func = eval(Symbol("win" * string(spectrum_type)))
	win_norm && (spec = win_norm_func(spec, window))

	info = LinSpecSetup(sr, freq_range, spectrum_type, win_norm)

	return LinSpec{typeof(stft)}(spec .* 2, freq, info)
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
Base.eltype(::LinSpec{T}) where T = T

get_data(s::LinSpec) = s.spec'
get_freq(s::LinSpec) = s.freq
get_info(s::LinSpec) = s.info
