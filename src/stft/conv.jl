# copiato da DSP.jl
# using FFTW

@inline function _zeropad!(
    padded::AbstractVector,
    u::AbstractVector
)
    data_region = CartesianIndices(u)
    data_dest = (1,)
    datasize = length(data_region)
    data_first_i = first(data_region)[1]
    dest_first_i = data_dest[1]
    copyto!(padded, dest_first_i, u, data_first_i, datasize)
    padded[1:dest_first_i-1] .= 0
    padded[dest_first_i+datasize:end] .= 0

    return padded
end

fast_fft_sizes = [2, 3, 5, 7]
nextfastfft(n) = nextprod(fast_fft_sizes, n)
nextfastfft(ns::Tuple) = nextfastfft.(ns)
nextfastfft(ns...) = nextfastfft(ns)

function conv(
    u::Union{AbstractVector{T},AbstractArray{T}},
    v::Union{AbstractVector{T},AbstractArray{T}},
    shape::Symbol
) where {T<:Real}

    su = size(u)
    sv = size(v)
    if prod(su) < prod(sv)
        u, v, su, sv = v, u, sv, su
    end

    outsize = su .+ sv .- 1
    out = similar(u, outsize)
    nffts = nextfastfft(outsize)

    padded = similar(u, nffts)
    _zeropad!(padded, u)
    p = plan_rfft(padded)
    uf = p * padded
    _zeropad!(padded, v)
    vf = p * padded
    uf .*= vf
    raw_out = irfft(uf, nffts[1])

    copyto!(
        out,
        CartesianIndices(out),
        raw_out,
        CartesianIndices(UnitRange.(1, outsize))
    )

    if (shape == :full)
        return out
    elseif (shape == :valid)
        valid_size = su .- sv .+ 1
        from = (outsize .- valid_size) .÷ 2 .+1
        to = from .+ valid_size .- 1

        return out[from[1]:to[1]]
    elseif (shape == :same)
        same_size = su
        from = Int.(round.((outsize .- same_size) ./ 2)) .+1
        to = from .+ same_size .- 1

        return out[from[1]:to[1]]
    else
        error("Unknown shape: '", shape, "'.")
    end
end # function conv

