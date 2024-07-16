# ---------------------------------------------------------------------------- #
#                            linear spectrogram                                #
# ---------------------------------------------------------------------------- #
# helper function to avoid log of zero
no_zero = x -> x == 0 ? floatmin(Float64) : x

function _get_spec(
    x::AbstractArray{T},
    # x_freq::StepRangeLen{Float64};
    x_freq;
    freq_range::Tuple{Int64, Int64},
    norm_func::Function,
    db_scale::Bool
) where T <: Union{Float64, Complex{Float64}}
    # trim to desired frequency range
    x_range = findall(freq_range[1] .<= x_freq .<= freq_range[2])
    lin_spec, lin_freq = x[x_range, :], x_freq[x_range]
    # lin_spec, lin_freq = x[:, x_range], x_freq[x_range]

    # apply normalization function
    lin_spec = norm_func.(lin_spec * 2)

    # scale to dB if requested
    if db_scale
        lin_spec = log10.(no_zero.(lin_spec))
    end

    return lin_spec, lin_freq
end

# ---------------------------------------------------------------------------- #
#                                  callings                                    #
# ---------------------------------------------------------------------------- #
function get_spectrogram!(
    rack::AudioRack,
    stft::Stft;
    freq_range::Tuple{Int64, Int64} = (0, floor(Int, rack.audio.sr / 2)),
    norm::Symbol = :power,
    db_scale::Bool = false
)
    # normalization functions
    norm_funcs = Dict(
        :none => x -> x,
        :power => x -> x / sum(rack.stft.win)^2,
        :magnitude => x -> x / sum(rack.stft.win),
    )

    # check if win_norm is valid
    @assert haskey(norm_funcs, norm) "Unknown win_norm: $win_norm."

    spec, freq = _get_spec(
        rack.stft.stft,
        rack.stft.freq;
        freq_range=freq_range,
        norm_func=norm_funcs[norm],
        db_scale=db_scale
    )

    rack.lin_spec = LinSpec(
        norm,
        db_scale,
        spec,
        freq
    )
end

function get_spectrogram!(
    rack::AudioRack,
    cwt::Cwt;
    freq_range::Tuple{Int64, Int64} = (0, floor(Int, rack.audio.sr / 2)),
    norm::Symbol = :power,
    db_scale::Bool = false
)
    # normalization functions
    norm_funcs = Dict(
        :none => x -> x,
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => x -> sqrt.(real.((x .* conj.(x))))
    )
    # check if win_norm is valid
    @assert haskey(norm_funcs, norm) "Unknown norm: $norm."

    spec, freq = _get_spec(
        rack.cwt.cwt,
        rack.cwt_fb.freq;
        freq_range=freq_range,
        norm_func=norm_funcs[norm],
        db_scale=db_scale
    )

    rack.cwt_spec = CwtSpec(
        norm,
        db_scale,
        spec,
        freq
    )
end