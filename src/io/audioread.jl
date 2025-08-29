function audioread(path::File{format"WAV"}; kwargs...)
    @show "wav"
end

function audioread(path::File{format"MP3"}; kwargs...)
    mpg123 = mpg123_new()
    mpg123_open(mpg123, filename(path))
    nframes = mpg123_length(mpg123)

    samplerate, nchannels, encoding = mpg123_getformat(mpg123)
    datatype = encoding_to_type(encoding)
    info = MP3Info(nframes, nchannels, samplerate, datatype)

    blocksize = mpg123_outblock(mpg123)
    encsize = sizeof(datatype)
    bufsize = div(blocksize, encsize * nchannels)
    # readbuf = Array{datatype, Int(nchannels)}(undef, nchannels, bufsize)
    source = MP3FileSource(filename(path), mpg123, info, bufsize)
    
    # buffer = try
    #     read(source)
    # finally
    #     close(source)
    # end
    # buffer

    read(source)
end

audioread(path::String; kwargs...) = audioread(query(path); kwargs...)