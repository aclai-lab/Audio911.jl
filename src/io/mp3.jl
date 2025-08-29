# ---------------------------------------------------------------------------- #
#                                mpg123 types                                  #
# ---------------------------------------------------------------------------- #
"""represents the C pointer mpg123_handle*. used by all mpg123 functions"""
const MPG123 = Ptr{Nothing}

const MPG123_OK         = 0
const MPG123_ERR        = -1
const MPG123_DONE       = -12
const MPG123_NEW_FORMAT = -11
const MPG123_NEED_MORE  = -10

# ---------------------------------------------------------------------------- #
#                               mp3 encodings                                  #
# ---------------------------------------------------------------------------- #
const MPG123_ENC_8           = 0x00f                                            # 0000 0000 0000 1111 Some 8 bit  integer encoding.
const MPG123_ENC_16          = 0x040                                            # 0000 0000 0100 0000 Some 16 bit integer encoding.
const MPG123_ENC_24          = 0x4000                                           # 0100 0000 0000 0000 Some 24 bit integer encoding.
const MPG123_ENC_32          = 0x100                                            # 0000 0001 0000 0000 Some 32 bit integer encoding.
const MPG123_ENC_SIGNED      = 0x080                                            # 0000 0000 1000 0000 Some signed integer encoding.
const MPG123_ENC_FLOAT       = 0xe00                                            # 0000 1110 0000 0000 Some float encoding.
const MPG123_ENC_SIGNED_16   = MPG123_ENC_16 | MPG123_ENC_SIGNED | 0x10         # 0000 0000 1101 0000 signed 16 bit
const MPG123_ENC_UNSIGNED_16 = MPG123_ENC_16 | 0x20                             # 0000 0000 0110 0000 unsigned 16 bit
const MPG123_ENC_UNSIGNED_8  = 0x01                                             # 0000 0000 0000 0001 unsigned 8 bit
const MPG123_ENC_SIGNED_8    = MPG123_ENC_SIGNED | 0x02                         # 0000 0000 1000 0010 signed 8 bit
const MPG123_ENC_ULAW_8      = 0x04                                             # 0000 0000 0000 0100 ulaw 8 bit
const MPG123_ENC_ALAW_8      = 0x08                                             # 0000 0000 0000 0100 alaw 8 bit
const MPG123_ENC_SIGNED_32   = MPG123_ENC_32 | MPG123_ENC_SIGNED | 0x1000       # 0001 0001 1000 0000 signed 32 bit
const MPG123_ENC_UNSIGNED_32 = MPG123_ENC_32 | 0x2000                           # 0010 0001 0000 0000 unsigned 32 bit
const MPG123_ENC_SIGNED_24   = MPG123_ENC_24 | MPG123_ENC_SIGNED | 0x1000       # 0101 0000 1000 0000 signed 24 bit
const MPG123_ENC_UNSIGNED_24 = MPG123_ENC_24 | 0x2000                           # 0110 0000 0000 0000 unsigned 24 bit
const MPG123_ENC_FLOAT_32    = 0x200                                            # 0000 0010 0000 0000 32bit float
const MPG123_ENC_FLOAT_64    = 0x400                                            # 0000 0100 0000 0000 64bit float

# ---------------------------------------------------------------------------- #
#                               pcm encodings                                  #
# ---------------------------------------------------------------------------- #
const PCM8Sample  = Fixed{Int8, 7}
const PCM16Sample = Fixed{Int16, 15}
const PCM20Sample = Fixed{Int32, 19}
const PCM24Sample = Fixed{Int32, 23}
const PCM32Sample = Fixed{Int32, 31}
const PCM64Sample = Fixed{Int64, 63}

# ---------------------------------------------------------------------------- #
#                                  MP3Info                                     #
# ---------------------------------------------------------------------------- #
struct MP3Info
    nframes    :: Int64
    nchannels  :: Int32
    samplerate :: Int64
    datatype   :: DataType

    function MP3Info(    
        nframes::Int64,
        nchannels::Int32,
        samplerate::Int64,
        datatype::DataType
    )
        new(nframes, nchannels, samplerate, datatype)
    end

    function MP3Info(buf::SampleBuf{T}) where {T}
        new(nframes(buf), nchannels(buf), samplerate(buf), T)
    end
end

# ---------------------------------------------------------------------------- #
#                                MP3FileSource                                 #
# ---------------------------------------------------------------------------- #
mutable struct MP3FileSource{T,N} <: AbstractSampleSource
    path    :: AbstractString
    mpg123  :: MPG123
    info    :: MP3Info
    pos     :: Int64
    readbuf :: Array{T,N}

    function MP3FileSource(path::AbstractString, mpg123::MPG123, info::MP3Info, bufsize::Integer)
        readbuf = Array{info.datatype,2}(undef, info.nchannels, bufsize)
        new{info.datatype,Int(info.nchannels)}(path, mpg123, info, Int64(0), readbuf)
    end
end

@inline nchannels(source::MP3FileSource)  = Int(source.info.nchannels)
@inline samplerate(source::MP3FileSource) = source.info.samplerate
@inline nframes(source::MP3FileSource)    = source.info.nframes
@inline Base.eltype(source::MP3FileSource{T}) where {T} = T

# ---------------------------------------------------------------------------- #
#                                mpg123 calls                                  #
# ---------------------------------------------------------------------------- #
"""create new mpg123 handle"""
function mpg123_new()
    err = Ref{Cint}(0)
    mpg123 = ccall((:mpg123_new, libmpg123), MPG123,
                   (Ptr{Cchar}, Ref{Cint}),
                   C_NULL, err)

    if err.x != MPG123_OK
        error("Could not create mpg123 handle: ", mpg123_plain_strerror(err.x))
    end

    mpg123
end

"""open an mp3 file at fiven path"""
function mpg123_open(mpg123::MPG123, path::AbstractString)
    err = ccall((:mpg123_open, libmpg123), Cint,
                (MPG123, Ptr{Cchar}),
                mpg123, path)

    if err != MPG123_OK
        mpg123_delete(mpg123)
        error("Could not open $path: ", mpg123_plain_strerror(err))
    end

    err
end

"""return the number of samples in the file"""
function mpg123_length(mpg123::MPG123)
    length = ccall((:mpg123_length, libmpg123), Int64, (MPG123,), mpg123)
    if length == MPG123_ERR
        error("Could not determine the frame length")
    end
    convert(Int64, length)
end

"""return birtate, number of channels and encoding of the mp3 file"""
function mpg123_getformat(mpg123::MPG123)
    bitrate = Ref{Clong}(0)
    nchannels = Ref{Cint}(0)
    encoding = Ref{Cint}(0)
    err = ccall((:mpg123_getformat, libmpg123), Cint,
                (MPG123, Ref{Clong}, Ref{Cint}, Ref{Cint}),
                mpg123, bitrate, nchannels, encoding)

    if err != MPG123_OK
        error("Could not read format: ", mpg123_plain_strerror(err))
    end

    bitrate.x, nchannels.x, encoding.x
end

"""return the appropriate block size for handling this mpg123 handle"""
function mpg123_outblock(mpg123::MPG123)
    ccall((:mpg123_outblock, libmpg123), Csize_t, (MPG123,), mpg123)
end

# ---------------------------------------------------------------------------- #
#                               mpg123 -> julia                                #
# ---------------------------------------------------------------------------- #
"""convert mpg123 encoding to julia datatype"""
function encoding_to_type(encoding)
    mapping = Dict{Integer, Type}(
       MPG123_ENC_SIGNED_16 => PCM16Sample,
       # TODO: support more
    )

    encoding in keys(mapping) || error("Unsupported encoding $encoding")
    mapping[encoding]
end