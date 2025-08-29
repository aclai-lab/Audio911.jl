module Audio911

using FileIO
# using SampledSignals: SampleBuf, SampleSource
# using WAV: wavread
# using DSP: resample

# ---------------------------------------------------------------------------- #
#                          FileIO Extension for MP3                            #
# ---------------------------------------------------------------------------- #
using mpg123_jll
# using LAME_jll
using FixedPointNumbers

include("io/signals.jl")

# MP3 files with ID3v1 (or no tags) and ID3v2 tags have different headers
mp3_headers = (UInt8[0xff, 0xfb], UInt8[0x49, 0x44, 0x33])
# add FileIO mp3 formats
try
    add_format(format"MP3", mp3_headers, [".mp3"], [:MP3])
catch
    @info "FileIO mp3 formats, aleady registered."
end

include("io/mp3.jl")

# ---------------------------------------------------------------------------- #
#                                   types                                      #
# ---------------------------------------------------------------------------- #
"""
    Maybe{T}

Type alias for `Union{T, Nothing}`.
"""
const Maybe{T} = Union{T, Nothing}

const MayInt = Maybe{Int}

include("io/audiofile.jl")
include("io/audioread.jl")
export audioread

end
