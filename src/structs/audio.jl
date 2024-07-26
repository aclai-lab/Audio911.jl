# ---------------------------------------------------------------------------- #
#                                    audio                                     #
# ---------------------------------------------------------------------------- #
struct Audio 
	data::AbstractVector{Float64}
	sr::Int64
end

# keyword constructor
function Audio(;
	data,
	sr
)
	audio = Audio(;
		data,
		sr
	)
	return audio
end

function Base.show(io::IO, audio::Audio)
	print(io, "Audio(data: $(length(audio.data)) samples, sr: $(audio.sr) Hz)")
end

function Base.display(audio::Audio)
    t = (0:length(audio.data)-1) / audio.sr
    p = plot(t, audio.data, 
             xlabel="Time (s)", 
             ylabel="Amplitude",
             title="Audio Waveform",
             legend=false)
    display(p)
end

function load_audio(;
    fname::AbstractString,
    sr::Int64 = 8000
)
    audio = Audio(
		py"load_audio"(fname, sr)...
	)
	return audio
end

function save_audio(;
	audio::Audio,
    fname::AbstractString
)
    py"save_audio"(fname, audio.data, audio.sr)
end