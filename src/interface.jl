"""
list of parameters used
# --- audio ------------------------------------------------------------------ #
sr,
norm = false
# --- stft ------------------------------------------------------------------- #
stft_length	= sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512,
win_type = (:hann, :periodic),
win_length = stft_length,
overlap_length = round(Int, stft_length / 2),
stft_norm = :power, # :none, :power, :magnitude, :pow2mag
"""

# ---------------------------------------------------------------------------- #
#                                 parameters                                   #
# ---------------------------------------------------------------------------- #
setup(; kwargs...) = NamedTuple(collect(kwargs))

f_prefix = "get_"
f_list = [:stft, :lin, :melfb, :mel, :mfcc, :deltas, :spec, :f0, :cwtfb, :cwt,]

# default_pipelines = Dict{Symbol,Vector{Symbol}}(
#     :stft -> [], 
#     :lin -> [:stft], 
#     :melfb -> [:stft], 
#     :mel -> [:stft, :melfb], 
#     :mfcc -> [:mel], 
#     :deltas -> [:stft, :melfb, :mel, :mfcc], 
#     :spec -> [:stft], 
#     :f0 -> [:stft], 
#     :cwtfb -> [], 
#     :cwt -> [:cwtfb]
# )

function audio911features(audio::Audio, feats::Union{Symbol, Tuple}, setup_params::NamedTuple)

end

function audio911features(
        file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, feats::Union{Symbol, Tuple}, setup_params::NamedTuple)
    if !isnothing(setup_params)
        sr, norm = get(setup_params, :sr, nothing), get(setup_params, :norm_audio, false)
    else
        sr, norm = nothing, false
    end
    audio911features(load_audio(; file = file, sr = sr, norm = norm), feats, setup_params)
end

audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, setup_params::NamedTuple, feats::Union{Symbol, Tuple}) = audio911features(file, feats, setup_params)

audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, feats::Union{Symbol, Tuple}) = audio911features(file, feats, NamedTuple())
audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, setup_params::NamedTuple) = audio911features(file, (), setup_params)

audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}) = audio911features(file, (), NamedTuple())