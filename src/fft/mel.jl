# using DSP

# include("../signalDataStructure.jl")
# include("fft.jl")

function hz2mel(
    hz::Vector{Float64},
    mel_style::Symbol=:htk # :htk, :slaney, :default_htk, :default_slaney
)
    if mel_style == :htk || mel_style == :default_htk
        mel = 2595 * log10.(1 .+ hz / 700)
    else # slaney
        linStep = 200 / 3
        logStep = log(6.4) / 27
        changePoint = 1000
        changePoint_mel = changePoint / linStep
        isLinearRegion = hz .< changePoint
        mel = deepcopy(hz)
        mel[isLinearRegion] .= hz[isLinearRegion] / linStep
        mel[.!isLinearRegion] .= changePoint_mel .+ log.(hz[.!isLinearRegion] / changePoint) / logStep
    end
    return mel
end # hz2mel

function mel2hz(
    mel::LinRange{Float64,Int64},
    mel_style::Symbol=:htk # :htk, :slaney, :default_htk, :default_slaney
)
    if mel_style == :htk || mel_style == :default_htk
        hz = 700 * (exp10.(mel / 2595) .- 1)
    else
        linStep = 200 / 3
        logStep = log(6.4) / 27
        changePoint = 1000
        changePoint_mel = changePoint / linStep
        isLinearRegion = mel .< changePoint_mel
        hz = [mel;]
        hz[isLinearRegion] .= hz[isLinearRegion] * linStep
        hz[.!isLinearRegion] .= changePoint * exp.(logStep * (mel[.!isLinearRegion] .- changePoint_mel))
    end
    return hz
end # mel2hz

function get_mel_norm_factor(spectrum_type::Symbol, fft_window::Vector{Float64})
    if spectrum_type == :power
        return 1 / (sum(fft_window)^2)
    else
        return 1 / sum(fft_window)
    end
end

### da generalizzare per 
### frequency_scale :mel, :bark, :erb
### filterbanl_design_domain :linear, :warped (da verificare se serve)
function designMelFilterBank(data::signal_data, setup::signal_setup)
    # set the design domain ### da implementare in futuro
    setup.filterbank_design_domain == :linear ? design_domain = :linear : design_domain = setup.frequency_scale

    # compute band edges
    # TODO da inserire il caso :erb e :bark
    melRange = hz2mel(setup.frequency_range, setup.mel_style)
    setup.band_edges = mel2hz(LinRange(melRange[1], melRange[end], setup.num_bands + 2), setup.mel_style)

    ### parte esclusiva per mel filterbank si passa a file designmelfilterbank.m
    # determine the number of bands
    num_edges = length(setup.band_edges)

    # determine the number of valid bands
    valid_num_edges = sum((setup.band_edges .- (setup.sr / 2)) .< sqrt(eps(Float64)))
    valid_num_bands = valid_num_edges - 2

    # preallocate the filter bank
    data.mel_filterbank = zeros(Float64, setup.fft_length, setup.num_bands)
    setup.mel_frequencies = setup.band_edges[2:end-1]

    # Set this flag to true if the number of FFT length is insufficient to
    # compute the specified number of mel bands
    FFTLengthTooSmall = false

    # if :hz 
    linFq = collect(0:setup.fft_length-1) / setup.fft_length * setup.sr

    # Determine inflection points
    @assert(valid_num_edges <= num_edges)
    p = zeros(Float64, valid_num_edges, 1)

    for edge_n in 1:valid_num_edges
        for index in eachindex(linFq)
            if linFq[index] > setup.band_edges[edge_n]
                p[edge_n] = index
                break
            end
        end
    end

    FqMod = linFq

    # Create triangular filters for each band
    bw = diff(setup.band_edges)

    for k in 1:Int(valid_num_bands)
        # Rising side of triangle
        for j in Int(p[k]):Int(p[k+1])-1
            data.mel_filterbank[j, k] = (FqMod[j] - setup.band_edges[k]) / bw[k]
        end
        # Falling side of triangle
        for j = Int(p[k+1]):Int(p[k+2])-1
            data.mel_filterbank[j, k] = (setup.band_edges[k+2] - FqMod[j]) / bw[k+1]
        end
        emptyRange1 = p[k] .> p[k+1] - 1
        emptyRange2 = p[k+1] .> p[k+2] - 1
        if (!FFTLengthTooSmall && (emptyRange1 || emptyRange2))
            FFTLengthTooSmall = true
        end
    end

    # mirror two sided
    range = get_onesided_fft_range(setup.fft_length)
    range = range[2:end]
    data.mel_filterbank[end:-1:end-length(range)+1, :] = data.mel_filterbank[range, :]

    data.mel_filterbank = data.mel_filterbank'

    # normalizzazione    
    BW = setup.band_edges[3:end] - setup.band_edges[1:end-2]

    if (setup.filterbank_normalization == :area)
        weight_per_band = sum(data.mel_filterbank, dims=2)
        if setup.frequency_scale != :erb
            weight_per_band = weight_per_band / 2
        end
    elseif (setup.filterbank_normalization == :bandwidth)
        weight_per_band = BW / 2
    else
        weight_per_band = ones(1, setup.num_bands)
    end

    for i = 1:setup.num_bands
        if (weight_per_band[i] != 0)
            data.mel_filterbank[i, :] = data.mel_filterbank[i, :] ./ weight_per_band[i]
        end
    end

    # get one side
    range = get_onesided_fft_range(setup.fft_length)
    data.mel_filterbank = data.mel_filterbank[:, range]
    # manca la parte relativa a :erb e :bark

    # setta fattore di normalizzazione
    if setup.window_norm
        win_norm_factor = get_mel_norm_factor(setup.spectrum_type, data.fft_window)
        data.mel_filterbank = data.mel_filterbank * win_norm_factor
    end
end # function designMelFilterBank

function createDCTmatrix(
    numCoeffs::Int64,
    numFilters::Int64,
    DT::DataType
)
    # create DCT matrix
    N = convert(DT, numCoeffs)
    K = numFilters
    matrix = zeros(DT, Int(N), numFilters)
    A = sqrt(1 / K)
    B = sqrt(2 / K)
    C = 2 * K
    piCCast = convert(DT, 2 * pi / C)

    matrix[1, :] .= A
    for k in 1:K
        for n in 2:Int(N)
            matrix[n, k] = B * cos(piCCast * (n - 1) * (k - 0.5))
        end
    end

    return matrix
end

function cepstralCoefficients(
    S::AbstractMatrix{T},
    numCoeffs::Int64,
    rectification::Symbol
) where {T<:AbstractFloat}
    DT = eltype(S)
    # Rectify
    if (rectification == :log)
        amin = floatmin(DT)
        S[S.==0] .= amin
        S = log10.(S)
    end
    # Reshape spectrogram to matrix for vectorized matrix multiplication
    L, M = size(S)
    N = 1
    S = reshape(S, L, M * N)
    # Design DCT matrix
    DCTmatrix = createDCTmatrix(numCoeffs, L, DT)
    # Apply DCT matrix
    coeffs = DCTmatrix * S

    return permutedims(coeffs, [2 1])

    # # Unpack matrix back to 3d array DA EVITARE già da prima
    # coeffs = reshape(coeffs,numCoeffs,M,N)
    # # Put time along the first dimension
    # coeffs = permutedims(coeffs,[2 1 3])
end # function cepstralCoefficients

function audioDelta(
    x::AbstractMatrix{T},
    window_length::Int64
) where {T<:AbstractFloat}

    DT = eltype(x)
    M = convert(DT, floor(window_length / 2))
    b = collect(M:-1:-M) ./ sum((1:M) .^ 2)

    return delta = filt(b, 1, x)
end

# function mel_spectrogram(
#     x::AbstractArray{T},
#     sr::Int64;
#     # default setup
#     # fft
#     window_type::Symbol=:hann,
#     window_length::Int=Int(round(0.03 * sr)),
#     overlap_length::Int=Int(round(0.02 * sr)),

#     # mel
#     num_bands::Int=32,
#     mel_style::Symbol=:slaney, # :htk, :slaney
#     frequency_range::Vector{Int64}=[0, Int(round(sr / 2))],
#     filterbank_normalization::Symbol=:bandwidth,
#     spectrum_type::Symbol=:power,

#     # filterbank_design_domain::Symbol=:linear, # settato, ma si usa?
#     # windowNormalization::Bool=true, # settato, ma si usa?
#     # oneSided::Bool=true # default, non viene parametrizzato
# ) where {T<:AbstractFloat}

#     # setup and data structures definition    
#     setup = signal_setup(
#         sr=sr,

#         # fft
#         window_type=window_type,
#         window_length=window_length,
#         overlap_length=overlap_length,

#         # linear spectrum
#         lin_frequency_range=[0.0, sr / 2],

#         # mel
#         num_bands=num_bands,
#         mel_style=mel_style,
#         frequency_range=Float64.(frequency_range),
#         filterbank_normalization=filterbank_normalization,
#         spectrum_type=spectrum_type,

#         # filterbank_design_domain=filterbank_design_domain, # settato, ma si usa?
#         # windowNormalization=windowNormalization, # settato, ma si usa?
#         # oneSided=oneSided # default, non viene parametrizzato
#     )

#     data = signal_data(
#         x=Float64.(x)
#     )

#     takeFFT(data, setup)
#     melSpectrogram(data, setup)
# end

function mel_spectrogram(
    data::signal_data,
    setup::signal_setup
)
    designMelFilterBank(data, setup)

    hop_length = setup.window_length - setup.overlap_length
    num_hops = Int(floor((size(data.x, 1) - setup.window_length) / hop_length) + 1)

    # apply filterbank
    if (setup.spectrum_type == :power)
        data.mel_spectrogram = reshape(data.mel_filterbank * data.fft, setup.num_bands, num_hops)
    else
        # magnitude, da fare
    end
    data.mel_spectrogram = data.mel_spectrogram'

end # melSpectrogram

function _mfcc(
    data::signal_data,
    setup::signal_setup
)

    # Calculate cepstral coefficients
    data.mfcc_coeffs = cepstralCoefficients(data.mel_spectrogram', setup.num_coeffs, setup.rectification)

    # place log energy
    if (setup.log_energy_pos == :append)
        data.mfcc_coeffs = hcat(data.mfcc_coeffs, data.log_energy)
    elseif (setup.log_energy_pos == :replace)
        # data.coeffs = [logE.',data.coeffs(:,2:end)];
        # da fare
    end

    # METTERE IL CASO CHE le delta vengono calcolate solo se necessario
    data.mfcc_delta = audioDelta(data.mfcc_coeffs, setup.delta_window_length)
    data.mfcc_deltadelta = audioDelta(data.mfcc_delta, setup.delta_window_length)

    # window_length = size(data.window, 1)
    # hop_length = window_length - setup.overlap_length
    # numHops = Int(floor((size(data.x, 1) - window_length) / hop_length) + 1)
    # bandUpLimit = collect(0:(numHops-1)) * hop_length .+ window_length

    # return data.mfcc_coeffs, data.mfcc_delta, data.mfcc_deltadelta, bandUpLimit
end

function mfcc(
    x::AbstractArray{T},
    sr::Int64;

    # default setup
    # fft
    window_type::Symbol=:hann,
    window_length::Int=Int(round(0.03 * sr)),
    overlap_length::Int=Int(round(0.02 * sr)),
    window_norm::Bool=true,

    # mel
    num_bands::Int=32,
    mel_style::Symbol=:slaney, # :htk, :slaney
    frequency_range::Vector{Int64}=[0, Int(round(sr / 2))],
    filterbank_normalization::Symbol=:bandwidth,
    spectrum_type::Symbol=:power,

    # mfcc
    num_coeffs::Int=13,
    rectification::Symbol=:log,
    log_energy_pos::Symbol=:append,
    delta_window_length::Int=9,

    # filterbank_design_domain::Symbol=:linear, # settato, ma si usa?
    # oneSided::Bool=true # default, non viene parametrizzato
) where {T<:AbstractFloat}

    # setup and data structures definition
    setup = signal_setup(
        sr=sr,

        # fft
        window_type=window_type,
        window_length=window_length,
        overlap_length=overlap_length,
        window_norm=window_norm,

        # linear spectrum
        lin_frequency_range=[0.0, sr / 2],

        # mel
        num_bands=num_bands,
        mel_style=mel_style,
        frequency_range=Float64.(frequency_range),
        filterbank_normalization=filterbank_normalization,
        spectrum_type=spectrum_type,

        # mfcc
        num_coeffs=num_coeffs,
        rectification=rectification,
        log_energy_pos=log_energy_pos,
        delta_window_length=delta_window_length,

        # filterbank_design_domain=filterbank_design_domain, # settato, ma si usa?
        # oneSided=oneSided # default, non viene parametrizzato
    )

    data = signal_data(
        x=Float64.(x)
    )

    takeFFT(data, setup)
    melSpectrogram(data, setup)
    _mfcc(data, setup)
end # mfcc