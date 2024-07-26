# ---------------------------------------------------------------------------- #
#                        continuous wavelets transform                         #
# ---------------------------------------------------------------------------- #
mutable struct Cwt
	# setup
    sr::Int64
	norm::Symbol
	db_scale::Bool
	# data
	spec::Union{Nothing, AbstractArray{Float64}}
	freq::Union{Nothing, AbstractVector{Float64}}
end

# keyword constructor
function Cwt(;
	sr,
	norm = :power, # :power, :magnitude, :pow2mag
    db_scale = false,
	spec = nothing,
	freq = nothing,
)
	cwt = Cwt(
        sr,
		norm,
        db_scale,
		spec,
		freq
	)
	return cwt
end

function _get_cwt(;
    x::AbstractVector{Float64}, 
    fbank::AbstractArray{Float64}, 
    freq::AbstractVector{Float64},
    cwt::Cwt
)
    # fourier transform of input
    # obtain the CWT in the Fourier domain
    # invert to obtain wavelet coefficients
    spec = ifft(fft(x)[1:size(fbank, 2)]' .* fbank, 2)

    # fft impleementation in matlab differs from julia
    # this code should be used for strong compatibility
    # but keep in mind that reversing the sign in imm part
    # will reverse the spectrogram
    # xposdft = fft(x)
    # cfsposdft = xposdft' .* fbank
    # cfsposdft_rev_imm = map(z -> real(z) - imag(z) * im, cfsposdft)
    # spec = ifft(cfsposdft_rev_imm, 2)

    # normalize
    norm_funcs = Dict(
        :power => x -> real.((x .* conj.(x))),
        :magnitude => x -> abs.(x),
        :pow2mag => x -> sqrt.(real.((x .* conj.(x))))
    )
    # check if spectrum_type is valid
    @assert haskey(norm_funcs, cwt.norm) "Unknown spectrum_type: $cwt.norm."

    cwt.spec = reverse(norm_funcs[cwt.norm](spec))
    cwt.freq = reverse(freq)

    return cwt
end

function Base.show(io::IO, cwt::Cwt)
    print(io, "Cwt(")
    print(io, "sr=$(cwt.sr), ")
    print(io, "norm=:$(cwt.norm), ")
    print(io, "db_scale=$(cwt.db_scale), ")
    if isnothing(cwt.spec)
        print(io, "spec=nothing, ")
    else
        print(io, "spec=$(size(cwt.spec)), ")
    end
    if isnothing(cwt.freq)
        print(io, "freq=nothing)")
    else
        print(io, "freq=$(length(cwt.freq)) frequencies)")
    end
end

function Base.display(cwt::Cwt)
    isnothing(cwt.spec) && error("Cwt has not been computed yet.")

    # apply log scaling to enhance contrast
    log_spec = log10.(cwt.spec .+ eps())
    
    # normalize to [0, 1] range
    normalized_spec = (log_spec .- minimum(log_spec)) ./ (maximum(log_spec) - minimum(log_spec))

    # define frequency ticks
    min_freq, max_freq = extrema(cwt.freq)
    freq_ticks = 10 .^ range(log10(min_freq), log10(max_freq), length=5)
    freq_ticks = round.(freq_ticks, digits=1)
    
    heatmap(1:size(cwt.spec, 2), cwt.freq, normalized_spec;
        ylabel="Frequency (Hz)",
        xlabel="Frame",
        title="CWT Spectrogram",
        color=:viridis,
        clim=(0, 1),
        yaxis=:log,
        yticks=(freq_ticks, string.(freq_ticks))
    )
end

function get_cwt(;
	audio::Audio,
    fbank::CwtFbank,
	kwargs...
)
	_get_cwt(x=audio.data, fbank=fbank.fbank, freq=fbank.freq, cwt=Cwt(;sr=audio.sr, kwargs...))
end