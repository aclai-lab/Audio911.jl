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
    file::Union{AbstractString, AbstractVector{Float64}},
    sr::Union{Nothing, Int64} = nothing,
	norm::Bool = false
)
	if file isa AbstractString
		audio = Audio(
			py"load_audio"(file, sr)...
		)
	elseif file isa AbstractVector{Float64} && sr isa Int64
		audio = Audio(file, sr)
	else
		throw(ArgumentError("Invalid arguments"))
	end

	# normalize audio
	if norm && length(audio.data) != 0
		audio.data = audio.data ./ maximum(abs.(audio.data))
	end

	return audio
end

function save_audio(;
	audio::Audio,
    file::AbstractString
)
    py"save_audio"(file, audio.data, audio.sr)
end