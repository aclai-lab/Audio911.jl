# ---------------------------------------------------------------------------- #
#                        continuous wavelets transform                         #
# ---------------------------------------------------------------------------- #
struct CwtSetup
	norm::Symbol
	dbscale::Bool
end

struct CwtData
    spec::AbstractArray{<:AbstractFloat}
	freq::AbstractVector{<:AbstractFloat}
end

struct Cwt
    sr::Int64
    setup::CwtSetup
    data::CwtData
end

function _get_cwt(;
    x::AbstractVector{Float64}, 
    sr::Int64,
    fbank::AbstractArray{Float64}, 
    freq::AbstractVector{Float64},
    norm::Symbol = :power, # :power, :magnitude, :pow2mag
    dbscale::Bool = false,
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
    @assert haskey(norm_funcs, norm) "Unknown spectrum_type: $norm."

    Cwt(sr, CwtSetup(norm, dbscale), CwtData(reverse(norm_funcs[norm](spec)), reverse(freq)))
end

function Base.show(io::IO, cwt::Cwt)
    print(io, "Cwt(")
    print(io, "sr=$(cwt.sr), ")
    print(io, "norm=:$(cwt.setup.norm), ")
    print(io, "db_scale=$(cwt.setup.dbscale), ")
end

function Base.display(cwt::Cwt)
    isnothing(cwt.data.spec) && error("Cwt has not been computed yet.")

    # apply log scaling to enhance contrast
    log_spec = log10.(cwt.data.spec .+ eps())
    
    # normalize to [0, 1] range
    normalized_spec = (log_spec .- minimum(log_spec)) ./ (maximum(log_spec) - minimum(log_spec))

    # define frequency ticks
    min_freq, max_freq = extrema(cwt.data.freq)
    freq_ticks = 10 .^ range(log10(min_freq), log10(max_freq), length=5)
    freq_ticks = round.(freq_ticks, digits=1)
    
    heatmap(1:size(cwt.data.spec, 2), cwt.data.freq, normalized_spec;
        ylabel="Frequency (Hz)",
        xlabel="Frame",
        title="CWT Spectrogram",
        color=:viridis,
        clim=(0, 1),
        yaxis=:log,
        yticks=(freq_ticks, string.(freq_ticks))
    )
end

function get_cwt(; source::Audio, fbank::CwtFb, kwargs...)
	_get_cwt(; x=source.data, sr=source.sr, fbank=fbank.data.fbank, freq=fbank.data.freq, kwargs...)
end