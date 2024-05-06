# reference https://www.audiolabs-erlangen.de/resources/MIR/chromatoolbox/
# reference https://www.ee.columbia.edu/~dpwe/resources/matlab/chroma-ansyn/#1

using Audio911
using DSP, LinearAlgebra, Statistics

function hz_to_octs(frequencies; tuning::Float64 = 0.0, bins_per_octave::Int = 12)
	A440 = 440.0 * 2.0^(tuning / bins_per_octave)
	octs = log2.(frequencies ./ (A440 / 16))
	return octs
end

function octs_to_hz(octs; tuning::Float64 = 0.0, bins_per_octave::Int = 12)
	A440 = 440.0 * 2.0^(tuning / bins_per_octave)
	return (A440 / 16) * (2.0 .^ octs)
end

# hz_to_octs = (frequencies, tuning=0.0, bins_per_octave=12) -> begin
#     A440 = 440.0 * 2.0 ^ (tuning / bins_per_octave)
#     log2.(frequencies ./ (A440 / 16))
# end

# octs_to_hz = (octs, tuning=0.0, bins_per_octave=12) -> begin
#     A440 = 440.0 * 2.0 ^ (tuning / bins_per_octave)
#     (A440 / 16) * (2.0 .^ octs)
# end

function A4_to_tuning(A4; bins_per_octave::Int = 12)
	tuning = bins_per_octave * (log2.(A4) - log2(440.0))
	return tuning
end

function tuning_to_A4(tuning; bins_per_octave::Int = 12)
	return 440.0 * 2.0^(tuning / bins_per_octave)
end

function design_chroma_filterbank(setup::AudioSetup, data::AudioData)
	#     sr::Float64,
	#     n_fft::Int,
	#     n_chroma::Int = 12,
	#     tuning::Float64 = 0.0,
	#     ctroct::Float64 = 5.0,
	#     octwidth::Union{Float64, Nothing} = 2.0,
	#     norm::Union{Float64, Nothing} = 2.0,
	#     base_c::Bool = true,
	# )
	# sr > setup.sr
	# n_fft > setup.fft_length
	# n_chroma > setup.bins_octave
	# tuning::Float64 = 0.0, deviazione da 440 da eliminare
	# ctroct::Float64 = 5.0, centroid frequency
	# octwidth::Union{Float64, Nothing} = 2.0, gaussian sd
	# norm::Union{Float64, Nothing} = 2.0, fattore di normalizzazione
	# base_c::Bool = true, parte da C, false parte da A

	wts = zeros(setup.bins_octave, setup.fft_length)

	# Get the FFT bins, not counting the DC component
	frequencies = range(0, sr, length = n_fft + 1)[2:end]

	frqbins = n_chroma * hz_to_octs(frequencies, tuning = tuning, bins_per_octave = n_chroma)

	# make up a value for the 0 Hz bin = 1.5 octaves below bin 1
	# (so chroma is 50% rotated from bin 1, and bin width is broad)
	frqbins = [frqbins[1] - 1.5 * n_chroma; frqbins]

	binwidthbins = max(frqbins[2:end] - frqbins[1:end-1], 1.0)

	D = reshape(repeat(frqbins, 1, n_chroma) .- repeat(0:n_chroma-1, length(frqbins)), n_chroma, length(frqbins))

	n_chroma2 = round(Int, n_chroma / 2)

	# Project into range -n_chroma/2.. n_chroma/2
	# add on fixed offset of 10*n_chroma to ensure all values passed to
	# rem are positive
	D = rem.(D .+ n_chroma2 .+ 10 * n_chroma, n_chroma) .- n_chroma2

	# Gaussian bumps - 2*D to make them narrower
	wts = exp.(-0.5 * (2 * D ./ repeat(binwidthbins, n_chroma, 1)) .^ 2)

	# normalize each column
	wts = normalize(wts, p = norm, dim = 2)

	# Maybe apply scaling for fft bins
	if !isnothing(octwidth)
		wts .*= repeat(exp.(-0.5 * (((frqbins / n_chroma .- ctroct) / octwidth) .^ 2)), n_chroma, 1)
	end

	if base_c
		wts = circshift(wts, (-3 * (n_chroma ÷ 12), 0))
	end

	# remove aliasing columns, copy to ensure row-contiguity
	return wts[:, 1:div(n_fft, 2)+1]
end

# @cache(level=10)
# def chroma(
#     *,
#     sr: float,
#     n_fft: int,
#     n_chroma: int = 12,
#     tuning: float = 0.0,
#     ctroct: float = 5.0,
#     octwidth: Union[float, None] = 2,
#     norm: Optional[float] = 2,
#     base_c: bool = True,
#     dtype: DTypeLike = np.float32,
# ) -> np.ndarray:

#     wts = np.zeros((n_chroma, n_fft))

#     # Get the FFT bins, not counting the DC component
#     frequencies = np.linspace(0, sr, n_fft, endpoint=False)[1:]

#     frqbins = n_chroma * hz_to_octs(
#         frequencies, tuning=tuning, bins_per_octave=n_chroma
#     )

#     # make up a value for the 0 Hz bin = 1.5 octaves below bin 1
#     # (so chroma is 50% rotated from bin 1, and bin width is broad)
#     frqbins = np.concatenate(([frqbins[0] - 1.5 * n_chroma], frqbins))

#     binwidthbins = np.concatenate((np.maximum(frqbins[1:] - frqbins[:-1], 1.0), [1]))

#     D = np.subtract.outer(frqbins, np.arange(0, n_chroma, dtype="d")).T

#     n_chroma2 = np.round(float(n_chroma) / 2)

#     # Project into range -n_chroma/2 .. n_chroma/2
#     # add on fixed offset of 10*n_chroma to ensure all values passed to
#     # rem are positive
#     D = np.remainder(D + n_chroma2 + 10 * n_chroma, n_chroma) - n_chroma2

#     # Gaussian bumps - 2*D to make them narrower
#     wts = np.exp(-0.5 * (2 * D / np.tile(binwidthbins, (n_chroma, 1))) ** 2)

#     # normalize each column
#     wts = util.normalize(wts, norm=norm, axis=0)

#     # Maybe apply scaling for fft bins
#     if octwidth is not None:
#         wts *= np.tile(
#             np.exp(-0.5 * (((frqbins / n_chroma - ctroct) / octwidth) ** 2)),
#             (n_chroma, 1),
#         )

#     if base_c:
#         wts = np.roll(wts, -3 * (n_chroma // 12), axis=0)

#     # remove aliasing columns, copy to ensure row-contiguity
#     return np.ascontiguousarray(wts[:, : int(1 + n_fft / 2)], dtype=dtype)

# def hz_to_octs(
#     frequencies: _ScalarOrSequence[_FloatLike_co],
#     *,
#     tuning: float = 0.0,
#     bins_per_octave: int = 12,
# ) -> Union[np.floating[Any], np.ndarray]:

#     A440 = 440.0 * 2.0 ** (tuning / bins_per_octave)

#     octs: np.ndarray = np.log2(np.asanyarray(frequencies) / (float(A440) / 16))
#     return octs

# def octs_to_hz(
#     octs: _ScalarOrSequence[_FloatLike_co],
#     *,
#     tuning: float = 0.0,
#     bins_per_octave: int = 12,
# ) -> Union[np.floating[Any], np.ndarray]:

#     A440 = 440.0 * 2.0 ** (tuning / bins_per_octave)

#     return (float(A440) / 16) * (2.0 ** np.asanyarray(octs))

# def A4_to_tuning(
#     A4: _ScalarOrSequence[_FloatLike_co], *, bins_per_octave: int = 12
# ) -> Union[np.floating[Any], np.ndarray]:

#     tuning: np.ndarray = bins_per_octave * (np.log2(np.asanyarray(A4)) - np.log2(440.0))
#     return tuning

# def tuning_to_A4(
#     tuning: _ScalarOrSequence[_FloatLike_co], *, bins_per_octave: int = 12
# ) -> Union[np.floating[Any], np.ndarray]:

#     return 440.0 * 2.0 ** (np.asanyarray(tuning) / bins_per_octave)




# using LinearAlgebra

# hz_to_octs = (frequencies, tuning=0.0, bins_per_octave=12) -> begin
#     A440 = 440.0 * 2.0 ^ (tuning / bins_per_octave)
#     map(f -> log2(f / (A440 / 16)), frequencies)
# end

# octs_to_hz = (octs, tuning=0.0, bins_per_octave=12) -> begin
#     A440 = 440.0 * 2.0 ^ (tuning / bins_per_octave)
#     map(o -> (A440 / 16) * (2.0 ^ o), octs)
# end

# A4_to_tuning = (A4, bins_per_octave=12) -> map(a -> bins_per_octave * (log2(a) - log2(440.0)), A4)

# tuning_to_A4 = (tuning, bins_per_octave=12) -> map(t -> 440.0 * 2.0 ^ (t / bins_per_octave), tuning)

# chroma = (; sr, n_fft, n_chroma=12, tuning=0.0, ctroct=5.0, octwidth=2, norm=2, base_c=true, dtype=Float32) -> begin
#     wts = zeros(n_chroma, n_fft)
#     frequencies = range(0, stop=sr, length=n_fft)[2:end]
#     frqbins = n_chroma .* hz_to_octs(frequencies, tuning=tuning, bins_per_octave=n_chroma)
#     frqbins = vcat([frqbins[1] - 1.5 * n_chroma], frqbins)
#     binwidthbins = vcat(max.(frqbins[2:end] .- frqbins[1:end-1], 1.0), [1])
#     D = [frqbins .- i for i in 0:n_chroma-1]'
#     n_chroma2 = round(n_chroma / 2)
#     D = map(d -> mod(d + n_chroma2 + 10 * n_chroma, n_chroma) - n_chroma2, D)
#     wts = map(d -> exp(-0.5 * ((2 * d / binwidthbins') ^ 2)), D)
#     wts = wts ./ sum(wts, dims=1)
#     if octwidth != nothing
#         wts .*= map(f -> exp(-0.5 * (((f / n_chroma - ctroct) / octwidth) ^ 2)), frqbins)'
#     end
#     if base_c
#         wts = circshift(wts, -3 * (div(n_chroma, 12)), dims=1)
#     end
#     return wts[:, 1:div(n_fft, 2)+1]
# end


#------------------------------------------------------------------------------#
#                                    debug                                     #
#------------------------------------------------------------------------------#
using Audio911
TESTPATH = joinpath(dirname(pathof(Audio911)), "..", "test")
TESTFILE = "common_voice_en_23616312.wav"
sr_src = 8000
x, sr = load_audio(joinpath(TESTPATH, TESTFILE), sr = sr_src)

setup = audio_setup(sr)
data = AudioData(Float64.(x))


