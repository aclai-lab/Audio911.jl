# ---------------------------------------------------------------------------- #
#                                    info                                      #
# ---------------------------------------------------------------------------- #
struct MelSpecSetup <: AbstractSetup
	sr       :: Int64
	win_norm :: Bool
end

# ---------------------------------------------------------------------------- #
#                          linear spectrogram struct                           #
# ---------------------------------------------------------------------------- #
struct MelSpec{F,T,B} <: AbstractSpectrogram
    spec  :: Matrix{T}
	fbank :: FBank{B}
    info  :: MelSpecSetup

    function MelSpec{F}(
        spec  :: Matrix{T},
		fbank :: FBank{B},
        info  :: MelSpecSetup
    ) where {F,T,B}
        new{F,T,B}(spec, fbank, info)
    end
end

# ---------------------------------------------------------------------------- #
#                                 get mel spec                                 #
# ---------------------------------------------------------------------------- #
function MelSpec(
	stft     :: Stft,
    fbank    :: AbstractFBank;
	win_norm :: Bool=false
)
    spectrum   = get_spectrum(stft)
	window     = get_window(stft)
	spec       = get_data(stft)
	filterbank = get_data(fbank)

	if win_norm
		spectrum == power     && (filterbank *= winpower(1, window))
		spectrum == magnitude && (filterbank *= winmagnitude(1, window))
	end

	info = MelSpecSetup(get_sr(stft), win_norm)

	return MelSpec{typeof(stft)}(filterbank * spec, fbank, info)
end

function MelSpec(
	stft     :: Stft;
    win_norm ::Bool=true,
    kwargs... # auditory filterbank kwargs
)
    fbank = auditory_fbank(get_sr(stft); sfreq=get_freq(stft), nfft=get_nfft(stft), kwargs...)
    MelSpec(stft, fbank; win_norm)
end

# ---------------------------------------------------------------------------- #
#                                    methods                                   #
# ---------------------------------------------------------------------------- #
Base.eltype(::MelSpec{T}) where T = T

get_data(m::MelSpec)  = m.spec'
get_freq(m::MelSpec)  = mem_freq(m.fbank)
get_setup(m::MelSpec) = m.info