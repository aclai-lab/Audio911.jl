using Test
using Audio911

test_files_dir()    = joinpath(dirname(@__FILE__), "test_files")
test_file(filename) = joinpath(test_files_dir(), filename)

wav_file = test_file("test.wav")
mp3_file = test_file("test.mp3")

audiofile = load(wav_file; mono=true, norm=false)

win = MovingWindow(size=512, step=256)
frames = AudioFrames(audiofile; win, type=hamming, periodic=true)

stft = Stft(frames; spectrum_type=magnitude, freq_range=(100,1000), win_norm=winmagnitude)

##############################################
function get_freq_range(
    frequency_range :: FreqRange,
    stft_size       :: Int64,
    sr              :: Int64
)
    # convert frequencies to bin indices
    bin_low = ceil(Int, get_low(frequency_range) * stft_size / sr + 1)
    bin_high = floor(Int, get_hi(frequency_range) * stft_size / sr + 1)

    bin_low = cld(get_low(frequency_range) * stft_size, sr) + 1
    bin_high = fld(get_hi(frequency_range) * stft_size, sr) + 1
    
    # create bit vector for the frequency range
    freq_mask = falses(stft_size >> 1 + 1)
    freq_mask[bin_low:bin_high] .= true

    return freq_mask
end

#------------------------------------------------------------------------------#
#                           spectrum normalizations                            #
#------------------------------------------------------------------------------#
function none() end
winpower(f, w)     = @. f / sum(w)^2 * 2
winmagnitude(f, w) = @. f / sum(w) * 2

# function LinSpec
# frequency_range :: Union{Tuple{Int64,Int64},FreqRange}=FreqRange(0, get_info(frames).sr>>1),
frequency_range=FreqRange(100,1000)
spec = stft.spec
freq = stft.freq
norm_factor = winmagnitude


frequency_range isa Tuple && (frequency_range = FreqRange(first(frequency_range), last(frequency_range)))

stft_size = get_info(stft).stft_size
sr = get_info(stft).sr

mask = get_freq_range(frequency_range, stft_size, sr)
norm_factor == none && (norm_factor(f, _) = @. f * 2)
spec_mask = norm_factor(spec[mask, :], )
freq_mask = freq[mask]

