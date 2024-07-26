# ---------------------------------------------------------------------------- #
#                                    stft                                      #
# ---------------------------------------------------------------------------- #
mutable struct Stft
	# setup
    sr::Int64
    x_length::Int64
	stft_length::Int64
	win_type::Tuple{Symbol, Symbol}
	win_length::Int64
	overlap_length::Int64
	norm::Symbol
	# data
	spec::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
    frames::Union{Nothing, AbstractArray{Float64}}
end

# keyword constructor
function Stft(;
	sr,
    x_length,
	stft_length	= sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512,
	win_type = (:hann, :periodic),
	win_length = stft_length,
	overlap_length = round(Int, stft_length / 2),
	norm = :power, # :none, :power, :magnitude, :pow2mag
	spec = nothing,
	freq = nothing,
    frames = nothing
)
	stft = Stft(
        sr,
        x_length,
		stft_length,
		win_type,
		win_length,
		overlap_length,
		norm,
		spec,
		freq,
        frames
	)
	return stft
end

function _get_stft(;
	x::AbstractVector{Float64},
	stft::Stft
)
	frames, _, wframes, _, _ = _get_frames(x, stft.win_type, stft.win_length, stft.overlap_length)

    if stft.win_length < stft.stft_length
        wframes = vcat(wframes, zeros(Float64, stft.stft_length - stft.win_length, size(wframes, 2)))
    elseif stft.win_length > stft.stft_length
        @error("FFT window size smaller than actual window size is highly discuraged.")
    end

    # take one side
    if mod(stft.stft_length, 2) == 0
        one_side = 1:Int(stft.stft_length / 2 + 1)   # even
    else
        one_side = 1:Int((stft.stft_length + 1) / 2)  # odd
    end

    # normalize
    norm_funcs = Dict(
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => x -> sqrt.(real.((x .* conj.(x))))
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, stft.norm) "Unknown spectrum_type: $stft.norm."

    stft.spec = norm_funcs[stft.norm](fft(wframes, (1,))[one_side, :])
    stft.freq = (stft.sr / stft.stft_length) * (one_side .- 1)
    stft.frames = frames

	return stft
end

function Base.show(io::IO, stft::Stft)
    print(io, "Stft(")
    print(io, "stft_length=$(stft.stft_length), ")
    print(io, "win_type=$(stft.win_type), ")
    print(io, "win_length=$(stft.win_length), ")
    print(io, "overlap_length=$(stft.overlap_length), ")
    print(io, "norm=:$(stft.norm), ")
    print(io, "spec=$(isnothing(stft.spec) ? "nothing" : size(stft.spec)), ")
    print(io, "freq=$(isnothing(stft.freq) ? "nothing" : length(stft.freq)))")
end

function Base.display(stft::Stft)
    isnothing(stft.spec) && error("STFT has not been computed yet.")
    
    heatmap(1:size(stft.spec, 2), stft.freq, stft.spec;
            ylabel="Frequency (Hz)",
            xlabel="Frame",
            title="STFT Spectrogram",
            color=:viridis,
			clim=(0, 1)
    )
end

function get_stft(;
	audio::Audio,
	kwargs...
)
	_get_stft(x=audio.data, stft=Stft(;sr=audio.sr, x_length=size(audio.data, 1), kwargs...))
end