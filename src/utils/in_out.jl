function load_audio(
    filename::AbstractString;
    sr::Int64=16000
)
    x, sr_def = py"load_audio"(filename, sr)
end

function save_audio(
    filename::AbstractString,
    x::AbstractVector{T};
    sr::Int64=16000
) where {T<:AbstractFloat}
    py"save_audio"(filename, x, sr)
end