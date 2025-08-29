# ---------------------------------------------------------------------------- #
#                               abstract types                                 #
# ---------------------------------------------------------------------------- #
abstract type AbstractAudioFile end

# ---------------------------------------------------------------------------- #
#                                 utilities                                    #
# ---------------------------------------------------------------------------- #
function convert2mono(audio::Matrix)::Vector{Abstractvector}
    return vec(mean(audio, dims=2))  # Average across channels
end

# ---------------------------------------------------------------------------- #
#                                 AudioFile                                    #
# ---------------------------------------------------------------------------- #
struct AudioFile <: AbstractAudioFile
    audio :: Vector{Float32}
    sr    :: Int
end

# ---------------------------------------------------------------------------- #
#                               AudioFile methods                                #
# ---------------------------------------------------------------------------- #
Base.length(w::AudioFile) = length(w.audio)

get_audio(w::AudioFile)::Vector{Float32} = w.audio
get_sr(w::AudioFile)::Int = w.sr
