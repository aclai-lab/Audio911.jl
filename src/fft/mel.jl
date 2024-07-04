# ---------------------------------------------------------------------------- #
#                               mel spectrogram                                #
# ---------------------------------------------------------------------------- #
function get_mel_spec!(
        setup::AudioSetup,
        data::AudioData
)
    filterbank = design_filterbank!(data, setup)

    hop_length = setup.stft.win_length - setup.stft.overlap_length
    num_hops = Int(floor((size(data.x, 1) - setup.stft.win_length) / hop_length) + 1)

    # apply filterbank
    # if (setup.spectrum_type == :power)
    data.mel_spectrogram = reshape(filterbank * data.stft.stft, setup.mel_bands, num_hops)
    # else
    #     #TODO
    #     error("magnitude not yet implemented.")
    # end

    data.mel_spectrogram = transpose(data.mel_spectrogram)
end # melSpectrogram

# ---------------------------------------------------------------------------- #
#                             mfcc related functions                           #
# ---------------------------------------------------------------------------- #
function create_DCT_matrix(
        mel_coeffs::Int64,
)
    # create DCT matrix
    matrix = zeros(Float64, mel_coeffs, mel_coeffs)
    s0 = sqrt(1 / mel_coeffs)
    s1 = sqrt(2 / mel_coeffs)
    piCCast = 2 * pi / (2 * mel_coeffs)

    matrix[1, :] .= s0
    for k in 1:mel_coeffs, n in 2:mel_coeffs
        matrix[n, k] = s1 * cos(piCCast * (n - 1) * (k - 0.5))
    end

    matrix
end

function audioDelta(
        x::AbstractMatrix{T},
        win_length::Int64,
        source::Symbol = :standard
) where {T <: AbstractFloat}

    # define window shape
    m = Int(floor(win_length / 2))
    b = collect(m:-1:(-m)) ./ sum((1:m) .^ 2)

    if source == :transposed
        filt(b, 1.0, x')'   #:audioflux setting
    else
        filt(b, 1.0, x)     #:matlab setting
    end
end

function get_mfcc!(
        setup::AudioSetup,
        data::AudioData
)
    # Rectify
    mel_spec = deepcopy(data.mel_spectrogram')

    if setup.normalization_type == :standard
        mel_spec[mel_spec .== 0] .= floatmin(Float64)
    elseif setup.normalization_type == :dithered
        mel_spec[mel_spec .< 1e-8] .= 1e-8
    else
        @warn("Unknown $setup.normalization_type normalization type, defaulting to standard.")
        mel_spec[mel_spec .== 0] .= floatmin(Float64)
    end

    # Design DCT matrix
    DCTmatrix = create_DCT_matrix(setup.mel_bands)

    # apply DCT matrix
    if (setup.rectification == :log)
        coeffs = DCTmatrix * log10.(mel_spec)
    elseif (setup.rectification == :cubic_root)
        # apply DCT matrix
        coeffs = DCTmatrix * mel_spec .^ (1 / 3)
    else
        @warn("Unknown $rectification DCT matrix rectification, defaulting to log.")
        coeffs = DCTmatrix * log10.(mel_spec)
    end

    # reduce to mfcc coefficients
    data.mfcc_coeffs = coeffs[1:(setup.mfcc_coeffs), :]'

    # log energy calc
    if setup.log_energy_source == :mfcc
        log_energy = sum(eachrow(mel_spec .^ 2)) / setup.mel_bands

        if setup.normalization_type == :standard
            log_energy[log_energy .== 0] .= floatmin(Float64)
        elseif setup.normalization_type == :dithered
            log_energy[log_energy .< 1e-8] .= 1e-8
        end

        data.log_energy = log.(log_energy)
    end

    if (setup.log_energy_pos == :append)
        data.mfcc_coeffs = hcat(data.mfcc_coeffs, data.log_energy)
    elseif (setup.log_energy_pos == :replace)
        data.mfcc_coeffs = hcat(data.log_energy, data.mfcc_coeffs[:, 2:end])
    end
end

function get_mfcc_deltas!(
        setup::AudioSetup,
        data::AudioData
)
    data.mfcc_delta = audioDelta(
        data.mfcc_coeffs, setup.delta_win_length, setup.delta_matrix)
    data.mfcc_deltadelta = audioDelta(
        data.mfcc_delta, setup.delta_win_length, setup.delta_matrix)
end
