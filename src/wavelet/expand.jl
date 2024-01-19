include("../fft/conv.jl")

if (!@isdefined(dwtAttr))
    struct dwtAttr{T<:Real}
        extMode::Symbol
        shift1D::Union{Int64,Float64}
        shift2D::AbstractVector{T}
    end
end # dwtAttr

function findactn(
    t::waveletPacketTree
) # find active nodes
    nodes = allnodes(t)
    type = :flag
    order = t.wavInfo.order
    spsch = t.tnodes.spsch
    nodes = depo2ind(order, nodes)
    tnrank = istnode(t, nodes)
    i_loc = locnumcn(order, nodes)
    act = spsch[Int.(i_loc)]
    i_Root = findfirst(x -> x == 0, nodes)

    if (!isempty(i_Root) && order > 0)
        act[i_Root] = true
    end

    for i in eachindex(tnrank)
        if (act[i] == 0)
            tnrank[i] = NaN
        end
    end

    return tnrank, nodes
end # function findactn

function nbnodes(
    order::Int64,
    depth::Int64,
)
    if (order == 0)
        return 0
    elseif (order == 1)
        return depth
    else
        return (order^(depth + 1) - 1) / (order - 1)
    end
end # function nbnodes

function wentropy(
    x::Union{AbstractVector{T},AbstractArray{T}},
    eName::Symbol,
    ep::Real
) where {T<:Real}

    if (eName == :shannon)
        x = x[x.!=0] .^ 2 # seleziona solo gli elementi != 0 e li eleva a 2
        ent = -sum(x .* log.(2^(-52) .+ x))
    elseif (eName == :threshold)

    elseif (eName == :sure)

    elseif (eName == :norm)

    elseif (eName == :energy || eName == :logenergy)

    end

    prec = 1.0e-10
    if (abs(ent) < prec)
        ent = 0
    end

    return ent
end # function wentropy

function defaninf(
    t::waveletPacketTree,
    nodes::Union{Real,AbstractMatrix{S}},
    x::Union{AbstractVecOrMat{T},AbstractArray{T}}
) where {T<:Real,S<:Real} #returns an array of numbers which contains information related to the nodes N

    nb = length(nodes)
    info = zeros(nb, 2)
    entname = t.entInfo.entName
    entpar = t.entInfo.entPar
    for k = 1:nb
        info[k, :] .= (wentropy(x[:, k], entname, entpar), NaN)
    end

    return info
end # function defaninf

function defaninf(t::waveletPacketTree, nodes::Union{Real,AbstractMatrix{S}}, x::Vector{Vector{T}}) where {T<:Real,S<:Real}
    defaninf(t, nodes, stack(x, dims=2))
end

function makeSym(
    x::Union{AbstractVector{T},AbstractArray{T}},
    lx::Int64,
    lf::Int64
) where {T<:Real}

    if (lx >= lf)
        xout = vcat(x[lf:-1:1], x[1:lx], x[lx:-1:lx-lf+1])
        xout = reshape(xout, 1, :)
    else

    end
end # function makeSym

function wextend(
    x::Union{AbstractVector{T},AbstractArray{T}},
    type::Symbol,
    mode::Symbol,
    lf::Int64
) where {T<:Real}

    if (type == :oneD)
        isOneD = 1
    end

    sz = (size(x, 1), size(x, 2))
    isROW = (sz[1] == 1)
    x = reshape(x, 1, :) # inverte le dimensioni
    sx = max(sz[1], sz[2])

    if (mode == :sym)
        x = makeSym(x, sx, lf)
    else

    end

    if (!isROW)
        x = vec(reshape(x, :, 1)) # inverte le dimensioni
    end
end # function

function dwt(
    x::Union{AbstractVector{T},AbstractArray{T}},
    lo_D::AbstractVector{S},
    hi_D::AbstractVector{S}
) where {T<:Real,S<:Real} # single-level discrete 1-D wavelet transform

    next = 3
    TFodd = mod(length(x), 2)
    dwt_attribute = dwtAttr(:sym, 0, [0, 0])
    dwtEXTM = dwt_attribute.extMode
    shift = dwt_attribute.shift1D

    # compute sizes and shape
    lx = length(x)
    lf = length(lo_D)
    f = 2 - shift

    if (dwtEXTM == :per)
        lenEXT = lf / 2
        l = 2 * ceil(lx / 2)
    else
        lenEXT = lf - 1
        l = lx + lf - 1
    end

    y = wextend(x, :oneD, dwtEXTM, lenEXT)
    # compute coefficients of approximation
    z = conv(y, lo_D, :valid)
    a = z[f:2:l]
    # compute coefficients of detail
    z = conv(y, hi_D, :valid)
    d = z[f:2:l]

    return a, d
end

function wsplit(
    t::waveletPacketTree,
    node::Real,
    x::Union{AbstractVector{T},AbstractArray{T}}
) where {T<:Real}  # split (decompose) the data of a terminal node
    order = t.wavInfo.order
    tnd = [Float64[] for i = 1:order] # quando in Matlab incontri cell (array of vectors)
    lo_D = t.wavInfo.lo_D
    hi_D = t.wavInfo.hi_D

    if (order == 2)
        tnd[1], tnd[2] = dwt(x, lo_D, hi_D)
    elseif (order == 4)
        tnd[1], tnd[2], tnd[3], tnd[4] = dwt2(x, lo_D, hi_D)
    end

    return tnd
end

function fmdtree(
    option::Symbol,
    t::waveletPacketTree,
    sizes::Vector{Any},
    tmpTMP::AbstractMatrix{T}
) where {T<:Real} #Field manager for DTREE object

    if (option == :setinit)

    elseif (option == :getinit)

    elseif (option == :an_del)

    elseif (option == :an_write)

    elseif (option == :an_read)

    elseif (option == :tn_beglensiz)

    elseif (option == :tn_write)
        # caso 4
        n = terminalNodes(sizes, tmpTMP)
        t.tN = n

    elseif (option == :tn_read)

    else
        error("Unknown option: '", option, "'.")
    end

    return t
end # function fmdtree

function expand(
    t::waveletPacketTree
)
    order = t.wavInfo.order
    depth = t.wavInfo.depth
    tnrank, nodes = findactn(t)
    n2dec = nodes[tnrank.==0]
    nbTOT = nbnodes(order, depth)
    x = (t.wavInfo.data)
    ndimX = ndims(x)
    last_COL_dim = ndimX + 1
    aninf = defaninf(t, 0, x)
    t.aN[1].ent, t.aN[1].ento = aninf
    data = [[] for i = 1:nbTOT]
    data[1] = x

    for j in eachindex(n2dec)
        node = n2dec[j]
        ind = node + 1
        x = Float64.(data[Int(ind)])
        tnval = wsplit(t, node, x)
        child = node * order .+ (1:order)'
        i_c = child .+ 1
        for k in 1:order
            data[Int(i_c[k])] = tnval[k]
            t.aN[Int(i_c[k])].size = (size(data[Int(i_c[k])], 1), size(data[Int(i_c[k])], 2))
        end
        data[Int(ind)] = []
        aninf = defaninf(t, child, tnval)
        for i in eachindex(i_c)
            t.aN[Int(i_c[i])].ent, t.aN[Int(i_c[i])].ento = aninf[i, :]
        end
    end

    ind_an = allnodes(t) .+ 1
    ind_tn = t.tnodes.tn .+ 1
    sizes = Vector(undef, length(ind_tn))

    for i in eachindex(ind_tn)
        sizes[i] = t.aN[Int.(ind_tn[i])].size[1], t.aN[Int.(ind_tn[i])].size[2]
    end

    lenTOT = sum(prod.(sizes))
    tmpTMP = zeros(1, lenTOT)
    iBEG = 1

    for k in eachindex(ind_tn)
        idx = ind_tn[k]
        iEND = iBEG + prod(sizes[k]) - 1
        tmpTMP[Int(iBEG):Int(iEND)] = data[Int(idx)]
        iBEG = iEND + 1
    end

    t = fmdtree(:tn_write, t, sizes, tmpTMP)

    return t
end # function expand