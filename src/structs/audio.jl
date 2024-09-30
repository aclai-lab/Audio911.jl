# ---------------------------------------------------------------------------- #
#                                    audio                                     #
# ---------------------------------------------------------------------------- #
struct Audio
    data::AbstractVector{Float64}
    sr::Int
end

const SUPPORTED_FORMATS = [".wav", ".mp3", ".ogg", ".flac"]

Base.show(io::IO, audio::Audio) = print(io, "Audio(length: $(length(audio.data)) samples, sr: $(audio.sr) Hz)")

function Plots.plot(audio::Audio; kwargs...)
    t = range(0, length = length(audio.data), step = 1/audio.sr)
    plot(t, audio.data; xlabel="Time (s)", ylabel="Amplitude", title="Audio Waveform", legend=false, kwargs...)
end

function load_audio(
        source::Union{AbstractString, AbstractVector{<:AbstractFloat}},
        sr::Union{Nothing, Int} = nothing;
        norm::Bool = false
)
    # check if stereo2mono is needed, TODO check
    validate_data = data -> begin
        if size(data, 2) > 1
            @warn "Multichannel audio detected: converting to mono."
            vec(mean(data, dims = 2))
		else
			data
		end
    end
    # scheck ample rate	
    validate_sr = sr -> begin
        if isnothing(sr)
            @warn "No sampling rate given, using default sampling rate of 8000Hz."
            8000
        else
            @assert sr > 0 && isinteger(sr) "Sample rate must be a positive integer."
            sr
        end
    end

    if source isa AbstractString
        # check file
        @assert isfile(source) "File does not exist."
        @assert any(endswith.(lowercase(source), SUPPORTED_FORMATS)) "Unsupported file format."

        data, sr = py"load_audio"(source, validate_sr(sr))
    else
        data = source
        if !isa(data, Vector{Float64})
            @warn "Converting audio to Float64."
            data = Float64.(data)
        end
    end

    # normalize audio
    if norm && !isempty(data)
        data ./= maximum(abs, data)
    end

    Audio(validate_data(data), validate_sr(sr))
end

function save_audio(audio::Audio, file::AbstractString)
    py"save_audio"(file, audio.data, audio.sr)
end
