function moving_mean(x::AbstractVector{Float64}, w::Int64)
    @assert isodd(w) "Window size must be odd"
    n = length(x)
    pad = w ÷ 2
    map(i -> median(@view x[max(1, i - pad):min(n, i + pad)]), 1:n)
end

function f_peaks(n::Vector{T}) where {T <: AbstractFloat}
    z7 = zeros(Float64, 7)
    z8 = zeros(Float64, 8)
    z3 = zeros(Float64, 3)

    nn = vcat(n, z7, n, z7, n, z8, n, z7, n, z7, n)

    n[end] = 0
    temp = repeat([z3; n; z3], 6, 1)

    b = all(reshape(nn .< temp, (Int(round(length(nn) / 6)), 6)), dims = 2)

    peaks_idx = []
    for i in eachindex(b)
        if b[i]
            push!(peaks_idx, i)
        end
    end
    peaks_idx = peaks_idx .- 3

    return peaks_idx
end

function get_threshs_from_feature(
        feature::Vector{Float64},
        bins::Int64,
        type::Symbol
)
    # get histogram
    hist_bins = round(Int, length(feature) / bins)
    # at leat 10 histogram
    hist_bins = max(10, hist_bins)

    m_feature = mean(feature)
    n_feature, edges_feature = get_histcounts(feature, nbins = hist_bins)

    # working with spectral spread
    if type == :specspread

        # delete the first indices of the histogram if the first bin contains 0s. We put them there when applying the threshold.
        if edges_feature[1] == 0
            n_feature = n_feature[2:end]
            edges_feature = edges_feature[2:end]
        end
        # set minimum values
        minval = m_feature / 2
    else
        minval = minimum(feature)
    end

    # find local maxima for energy
    peaks_idx = f_peaks(Float64.(n_feature))

    # working with energy
    if type == :energy
        # check to make sure that no peaks are beyond halfway through the histogram.
        # high maxima values later in the histogram makes very high threshold values.
        if length(peaks_idx) >= 2 && all(maximum(n_feature[peaks_idx[1:2]]) > m_feature)
            # Remove all maxima that are far along
            if edges_feature[peaks_idx[2]] > m_feature
                peaks_idx = peaks_idx[1]
            elseif edges_feature[peaks_idx[1]] > m_feature
                deleteat!(peaks_idx, 1)
            end
        end
    end

    if length(peaks_idx) == 0
        M1 = m_feature / 2
        M2 = minval
    elseif length(peaks_idx) == 1
        eF0 = vcat(collect(edges_feature), 0)
        M1 = 0.5 * (vcat(0, collect(edges_feature)) .- eF0) + eF0
        M1 = M1[peaks_idx .+ 1]
        M2 = minval
    else
        eF0 = vcat(collect(edges_feature), 0)
        AA = 0.5 * (vcat(0, collect(edges_feature)) .- eF0) + eF0
        M2 = AA[peaks_idx[1] + 1]
        M1 = AA[peaks_idx[2] + 1]
    end

    return M1, M2
end

function debuffer_frame_overlap(
        speech_mask::BitVector, window_length::Int64, overlap_length::Int64)
    hop_length = window_length - overlap_length

    n_shared_frames = floor(Int, window_length / hop_length)

    nearest_nv = DSP.filt(ones(n_shared_frames), [speech_mask; zeros(n_shared_frames - 1)])

    begin_thr = (1:(n_shared_frames - 1)) ./ 2
    end_thr = reverse(begin_thr)
    mid_thr = ones(length(nearest_nv) - 2 * (n_shared_frames - 1))

    thresh = [begin_thr; mid_thr; end_thr]

    out = nearest_nv .>= thresh

    return out, hop_length
end

function speech_detector(
        x::AbstractVector{Float64},
        sr::Int64;
        window_type::Tuple{Symbol, Symbol} = (:hann, :periodic),
        window_length::Int64 = round(Int, 0.03 * sr),
        overlap_length::Int64 = 0,
        thresholds::Tuple{Float64, Float64} = (-Inf, -Inf),
        merge_distance::Int64 = window_length * 5
)
    # ---------------------------------------------------------------------------------- #
    #                                     parameters                                     #
    # ---------------------------------------------------------------------------------- #
    # options = SdOptions(threshold, merge_distance)

    weight = 5 # weight for finding local maxima
    bins = 15 # number of bins for histograms
    spread_threshold = 0.05 # threshold used for not counting spectral spread when under this energy
    lower_spread_threshold_factor = 0.8 # factor to lower the spectral spread threshold.
    smoothing_filter_length = 5 # after getting spectral spread and energy data, 
    # these features are smoothed by filters of this length

    # ---------------------------------------------------------------------------------- #
    #      step 1: extract short-term spectral spread and energy from whole signal       #
    # ---------------------------------------------------------------------------------- #

    # normalize
    x = normalize_audio(x)

    # get stft
    frames, window = _get_frames(
        x,
        window_type = window_type,
        window_length = window_length,
        overlap_length = overlap_length
    )
    s, sfreq = _get_stft(
        frames .* window,
        sr,
        fft_length = 2 * window_length,
        window = window,
        frequency_range = (0, floor(Int, sr / 2)),
        spectrum_type = :magnitude
    )

    # determine short term energy
    energy = vec(window' .^ 2 * frames .^ 2)

    # filter the short term energy twice
    filtered_energy = moving_mean(
        moving_mean(energy, smoothing_filter_length), smoothing_filter_length)

    # get spectral spread and normalize
    spread = _get_spec_spread(s, sfreq) / (sr / 2)

    # set spectral spread value to 0 for frames with low energy
    spread[energy .< spread_threshold] .= 0

    # filter spectral spread twice
    filtered_spread = moving_mean(
        moving_mean(spread, smoothing_filter_length), smoothing_filter_length)

    #----------------------------------------------------------------------------------#
    #                        step 2: determine thresholds                              #
    #----------------------------------------------------------------------------------#
    # if both energy and spectral spread thresholds are defined
    if thresholds != (-Inf, -Inf)
        energy_thresh = thresholds[1]
        spread_thresh = thresholds[2]
    else
        # calculate energy and spectral spread local bin maxima
        e_m1, e_m2 = get_threshs_from_feature(filtered_energy, bins, :energy)
        s_m1, s_m2 = get_threshs_from_feature(filtered_spread, bins, :specspread)

        # calculate energy and spectral spread thresholds
        ww = 1 / (weight + 1)
        sspread_thresh = ww * (weight * s_m2 + s_m1[1]) * lower_spread_threshold_factor
        energy_Thresh = ww * (weight * e_m2 + e_m1[1])
    end
    #----------------------------------------------------------------------------------#
    #                       step 3: apply threshold criterion                          #
    #----------------------------------------------------------------------------------#
    speech_mask = (filtered_spread .> sspread_thresh) .& (filtered_energy .> energy_Thresh)

    #----------------------------------------------------------------------------------#
    #       step 4: Merge frames if they are close as described by MergeDistance       #
    #----------------------------------------------------------------------------------#
    # De-buffer for Frame Overlap
    if overlap_length > 0
        unbuff_out, window_length = debuffer_frame_overlap(
            speech_mask, window_length, overlap_length)
    else
        unbuff_out = speech_mask
    end

    # change frames into data points
    a = repeat(unbuff_out', outer = [window_length, 1])
    unbuff_out_mask = [a[:]; falses(length(x) - length(a), 1)]
    difference = diff([unbuff_out_mask; false], dims = 1)

    # find all changes from speech to silence; return index before change
    idx_m1 = findall(difference .== -1)
    idx_m1 = getindex.(idx_m1, 1)

    # find all changes from silence to speech; return index before change
    if unbuff_out[1] == 0
        idx_p1 = findall(difference .== 1)
        idx_p1 = getindex.(idx_p1, 1)
    else
        idx_p1 = findall(difference .== 1)
        idx_p1 = vcat(1, getindex.(idx_p1, 1))
    end

    # find gaps less than merge distance
    if length(idx_p1) > 1
        testmask = idx_p1[2:end] .- idx_m1[1:(length(idx_p1) - 1)] .<= merge_distance
    else
        testmask = falses(0, 1)
    end

    if isempty(idx_p1) || isempty(idx_m1)
        # no speech is found
        outidx = []
    else
        # arrange output
        idx_p2 = idx_p1[2:end, :]
        idx_m2 = idx_m1[1:(length(idx_p1) - 1), :]
        amask = .!testmask
        outidx = reshape([idx_p1[1]; idx_p2[amask]; idx_m2[amask]; idx_m1[end]], :, 2)
    end

    y = Float64[]
    for i in eachrow(outidx)
        y = [y; x[i[1]:i[2]]]
    end

    return y, outidx
end

function speech_detector(x::AbstractVector{<:AbstractFloat}, sr::Int64; kwargs...)
    speech_detector(Float64.(x), sr; kwargs...)
end

# references
# https://github.com/linan2/Voice-activity-detection-VAD-paper-and-code