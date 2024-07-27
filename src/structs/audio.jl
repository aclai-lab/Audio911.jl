# ---------------------------------------------------------------------------- #
#                                    audio                                     #
# ---------------------------------------------------------------------------- #
mutable struct Audio 
	data::AbstractVector{Float64}
	sr::Int64
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
    fname::Union{AbstractString, AbstractVector{Float64}},
    sr::Union{Nothing, Int64} = nothing,
	norm::Bool = false
)
	if fname isa AbstractString
		audio = Audio(
			py"load_audio"(fname, sr)...
		)
	elseif fname isa AbstractVector{Float64} && sr isa Int64
		audio = Audio(fname, sr)
	else
		throw(ArgumentError("Invalid arguments"))
	end

	# normalize audio
	if norm && length(audio.data) != 0
		audio.data ./ maximum(abs.(audio.data))
	end

	return audio
end

function save_audio(;
	audio::Audio,
    fname::AbstractString
)
    py"save_audio"(fname, audio.data, audio.sr)
end