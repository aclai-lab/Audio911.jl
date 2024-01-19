# lavora solo su segnali mono, quindi bisogna caricare i wav sempre con librosa > mono=true
using Statistics

include("../windowing/windowing.jl")

function signal_to_rms(
    y::Union{AbstractVector{T},AbstractArray{T}}
) where {T<:Real}
    # calcola power e rms
    return sqrt.(mean(abs.(y) .^ 2, dims=2))
end

function signal_to_db(
    y::Union{AbstractVector{T},AbstractArray{T}}
) where {T<:Real}
    magnitude = abs.(y)
    power = magnitude .^ 2

    ref_value = maximum(magnitude)^2
    amin = 1e-5^2

    # power to db
    db_log = 10.0 * log10.(max.(maximum.(power), amin))
    db_log .-= 10.0 * log10(max(amin, ref_value))
    return db_log
end

function trimAudio(
    x::Union{AbstractVector{T},AbstractArray{T}},
    fftLength::Int64=1024,
    threshold::Real=60
) where {T<:Real}

    # crea la matrice delle finestre: le righe sono le finestre.
    y = windowing(x, fftLength, :rect, :symmetric, false)

    # converte l'audio in rms e poi in db
    rms = signal_to_rms(y)
    db = signal_to_db(rms)
    framesValidated = db .> -threshold

    signal = zeros(eltype(x), 0)
    silence = zeros(eltype(x), 0)
    flag = :start
    for i in eachindex(framesValidated)
        # caso 1: partenza segnale
        if (framesValidated[i] && (flag == :silence || flag == :start))
            fade_y = fade(y[i, :], fftLength, :in)
            signal = vcat(signal, y[i, :])

            if (flag == :silence)
                silence = fade(silence, fftLength, :out)
            end

            flag = :signal
        # caso 2: continuazione segnale
        elseif (framesValidated[i] && flag == :signal)
            signal = vcat(signal, y[i, :])
        #caso 3: partenza silenzio
        elseif (!framesValidated[i] && (flag == :signal || flag == :start))
            fade_y = fade(y[i, :], fftLength, :in)
            silence = vcat(silence, y[i, :])

            if (flag == :signal)
                signal = fade(signal, fftLength, :out)
            end

            flag = :silence
            #caso 4: continuazione silenzio
        else
            silence = vcat(silence, y[i, :])
        end

        #chiudo con un fade out
        if (flag == :silence)
            silence = fade(silence, fftLength, :out)
        else
            signal = fade(signal, fftLength, :out)
        end
    end

    return signal, silence
end

##MAIN
# if (!@isdefined(mono)) # carica solo se non è definita mono, serve giusto per debug
#     # cosi ogni volta che lancio non ricarica tutto
#     using PyCall
#     librosa = pyimport("librosa")
#     sr_src = 8000
#     mono, sr = librosa.load("/home/riccardopasini/Documents/Julia/test_mono.wav", sr=sr_src, mono=true)
# end

# x, y = trimAudio(mono, 1024, 15)

# soundfile = pyimport("soundfile")
# ## Save audiofile
# soundfile.write(
#     "signal.wav",
#     x,
#     samplerate=sr,
#     subtype="PCM_16")
# soundfile.write(
#     "silence.wav",
#     y,
#     samplerate=sr,
#     subtype="PCM_16")


