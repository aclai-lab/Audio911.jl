"""
list of parameters used
# --- audio ------------------------------------------------------------------ #
sr,
norm = false
# --- stft ------------------------------------------------------------------- #
nfft	= sr !== nothing ? (sr <= 8000 ? 256 : 512) : 512,
win_type = (:hann, :periodic),
win_length = nfft,
overlap_length = round(Int, nfft / 2),
stft_norm = :power, # :none, :power, :magnitude, :pow2mag
"""

# ---------------------------------------------------------------------------- #
#                                 parameters                                   #
# ---------------------------------------------------------------------------- #
setup(; kwargs...) = NamedTuple(collect(kwargs))

f_prefix = "get_"
f_list = [:stft, :lin, :melfb, :mel, :mfcc, :deltas, :spec, :f0, :cwtfb, :cwt,]

# default_pipelines = Dict{Symbol,Vector{Symbol}}(
#     :stft -> [], 
#     :lin -> [:stft], 
#     :melfb -> [:stft], 
#     :mel -> [:stft, :melfb], 
#     :mfcc -> [:mel], 
#     :deltas -> [:stft, :melfb, :mel, :mfcc], 
#     :spec -> [:stft], 
#     :f0 -> [:stft], 
#     :cwtfb -> [], 
#     :cwt -> [:cwtfb]
# )

function audio911features(audio::Audio, feats::Union{Symbol, Tuple}, setup_params::NamedTuple)

end

function audio911features(
        file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, feats::Union{Symbol, Tuple}, setup_params::NamedTuple)
    if !isnothing(setup_params)
        sr, norm = get(setup_params, :sr, nothing), get(setup_params, :norm_audio, false)
    else
        sr, norm = nothing, false
    end
    audio911features(load_audio(; file = file, sr = sr, norm = norm), feats, setup_params)
end

audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, setup_params::NamedTuple, feats::Union{Symbol, Tuple}) = audio911features(file, feats, setup_params)

audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, feats::Union{Symbol, Tuple}) = audio911features(file, feats, NamedTuple())
audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}, setup_params::NamedTuple) = audio911features(file, (), setup_params)

audio911features(file::Union{AbstractString, AbstractVector{<:AbstractFloat}}) = audio911features(file, (), NamedTuple())

####################################################################################################################
####################################################################################################################
####################################################################################################################

function afe(x::Union{String, AbstractVector{Float64}}; kwargs...)
    p_mapping = Dict(:sr => :sr, :audio_norm => :norm)
    p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
    audio = load_audio(; source=x, p_kwargs...)

    featset = get(kwargs, :featset, ()) |> x -> x isa Symbol ? (x,) : x

    p_mapping = Dict(:nfft => :nfft, :fft_wintype => :wintype, :fft_winlength => :winlength, :fft_overlaplength => :overlaplength, :fft_norm => :norm)
    p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
    stftspec = get_stft(; source=audio, p_kwargs...)

    :lin in featset && begin
        p_mapping = Dict(:lin_freqrange => :freqrange, :lin_dbscale => :dbscale)
        p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
        linspec = get_linspec(; source=stftspec, p_kwargs...)
    end

    p_mapping = Dict(:mel_nbands => :nbands, :mel_scale => :scale, :mel_norm => :norm, :mel_freqrange => :freqrange, :mel_semitonerange => :semitonerange)
    p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
    melfb = get_melfb(; source=stftspec, p_kwargs...)

    :get_only_freqs in featset && return melfb.data.freq

    p_mapping = Dict(:mel_dbscale => :dbscale)
    p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
    melspec = get_melspec(; source=stftspec, fbank=melfb, p_kwargs...)

    :mfcc in featset && begin
        p_mapping = Dict(:mfcc_ncoeffs => :ncoeffs, :mfcc_rect => :rect, :mfcc_dither => :dither)
        p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
        mfcc = get_mfcc(; source=melspec, p_kwargs...)

        :deltas in featset && begin
            p_mapping = Dict(:delta_length => :dlength, :delta_transpose => :transpose)
            p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
            deltas = get_deltas(; source=mfcc, p_kwargs...)
        end
    end

    :f0 in featset && begin
        p_mapping = Dict(:f0_method => :method, :f0_freqrange => :freqrange, :f0_mflength => :mflength)
        p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
        f0 = get_f0(; source=stftspec, p_kwargs...)
    end

    p_mapping = Dict(:f0_method => :method, :f0_freqrange => :freqrange, :f0_mflength => :mflength)
    p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
    spect = get_spectrals(; source=stftspec, p_kwargs...)

    :cwt in featset && begin
        p_mapping = Dict(:cwt_wavelet => :wavelet, :cwt_morseparams => :morseparams, :cwt_vpo => :vpo, :cwt_freqrange => :freqrange)
        p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
        cwtfb = get_cwtfb(; source=audio, p_kwargs...)

        p_mapping = Dict(:cwt_norm => :norm, :cwt_dbscale => :dbscale)
        p_kwargs = Dict(p_mapping[k] => v for (k, v) in kwargs if k in keys(p_mapping))
        cwtspec = get_cwt(; source=audio, fbank=cwtfb, p_kwargs...)
    end

    hcat(
        filter(!isnothing, [
            :lin in featset ? linspec.data.spec' : nothing,
            melspec.data.spec',
            :mfcc in featset ? mfcc.data.spec' : nothing,
            :deltas in featset ? hcat(deltas.data.dspec', deltas.data.ddspec') : nothing,
            :f0 in featset ? f0.data.f0 : nothing,
            spect.data.centroid,
            spect.data.crest,
            spect.data.entropy,
            spect.data.flatness,
            spect.data.flux,
            spect.data.kurtosis,
            spect.data.rolloff,
            spect.data.skewness,
            spect.data.decrease,
            spect.data.slope,
            spect.data.spread
        ])...
    )  
end