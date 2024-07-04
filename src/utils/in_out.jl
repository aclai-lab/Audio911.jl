function load_audio(
    filename::AbstractString,
    sr::Int64=16000
)
    Audio(py"load_audio"(filename, sr)...)
end

function save_audio(
    filename::AbstractString,
    x::AbstractVector{<:AbstractFloat},
    sr::Int64=16000
)
    py"save_audio"(filename, x, sr)
end