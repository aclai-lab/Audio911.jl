function calc_window(
    half::Real,
    fftLength::Int64,
    winType::Symbol
)
    x = (0:half-1)' / (fftLength - 1)

    if (winType == :hann)
        w = 0.5 .- 0.5 * cos.(2 * pi * x)
    elseif (winType == :hamming)
        w = 0.54 .- 0.46 * cos.(2 * pi * x)
    elseif (winType == :blackman)
        w = 0.42 .- 0.5 * cos.(2 * pi * x) .+ 0.08 * cos.(4 * pi * x)
        if (!isempty(w))
            w[1] = 0
        end
    elseif (winType == :flattopwin)
        a0 = 0.21557895
        a1 = 0.41663158
        a2 = 0.277263158
        a3 = 0.083578947
        a4 = 0.006947368
        w = a0 .- a1 * cos.(2 * pi * x) .+ a2 * cos.(4 * pi * x) .- a3 * cos.(6 * pi * x) .+ a4 * cos.(8 * pi * x)
    end

    return w
end # function calc_window

function sym_window(
    fftLength::Int64,
    winType::Symbol,
    startIdx::Int
)
    if iseven(fftLength)
        half = fftLength / 2
        w = vec(calc_window(half, fftLength, winType)) # calcola la finestra dal min al max
        w = vcat(w, reverse(w[startIdx:end])) # completa la finestra appendendo la finestra specchiata
    else
        half = (fftLength + 1) / 2
        w = vec(calc_window(half, fftLength, winType)) # calcola la finestra dal min al max
        w = vcat(w, reverse(w[startIdx:end-1])) # completa la finestra appendendo la finestra specchiata
    end
    return w, half
end # function sym_window

function gencoswin(
    winType::Symbol,
    fftLength::Int64,
    samplingFlag::Symbol
)

    if (samplingFlag == :periodic)
        w, hL = sym_window(fftLength + 1, winType, 2)
    elseif (samplingFlag == :symmetric)
        w, hL = sym_window(fftLength, winType, 1)
    end
    return w, Int(round(hL))
end # function gencoswin

function rectwin(
    fftLength::Int64
)
        return ones(fftLength), fftLength
end