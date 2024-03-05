include("wpdec.jl")

"""
mutable struct allNodes{T<:Real}
    ind::Int64 # index
    size::Tuple{Int64,Int64} # size of data
    ent::T # entropy
    ento::T # optimal entropy
end

mutable struct terminalNodes{T<:Real}
    size::AbstractVector{Any}
    data::AbstractMatrix{T}
end

struct wTnodes
    tn::Vector{Float64} # Array of terminal nodes of tree
    nbtn::Int64
    spsch::Vector{Int64}
end

struct wavInfo{T<:Real,S<:Real}
    data::Union{AbstractVector{T},AbstractArray{T}}
    dataSize::Tuple{Int64,Int64} # size of initial data
    order::Int64
    depth::Int64
    wavName::String # wavelet name
    lo_D::AbstractVector{S} # low decomposition filter
    hi_D::AbstractVector{S} # high decomposition filter
    lo_R::AbstractVector{S} # low reconstruction filter
    hi_R::AbstractVector{S} # high reconstruction filter 
end

struct entInfo{T<:Real}
    entName::Symbol # entropy name
    entPar::T # entropy parameter
end

mutable struct waveletPacketTree
    aN::Vector{allNodes}
    tN::terminalNodes
    tnodes::wTnodes
    wavInfo::wavInfo
    entInfo::entInfo
end
"""

function ind2depo(
    order::Int64,
    node::AbstractVector{T}
) where {T<:AbstractFloat}

    r = size(node, 1)
    c = size(node, 2)

    if (c <= 1)
        d = zeros(r, 1)
        p = zeros(r, 1)
        K = node .> 0

        if any(K)
            if (order == 1)
                d[k] = node[k]
            else
                d[K] = floor.(log.((order - 1) * node[K] .+ 1) / log(order))
                p[K] = node[K] .- (order .^ d[K] .- 1) / (order - 1)
            end
        end

        d[node.<0] .= -1 # funziona

    elseif (c == 2)
        # ind2depo linea 45
    end

    return d, p
end

function frqord(
    node::AbstractVector{T}
) where {T<:AbstractFloat}

    order = 2
    depths, pos_nat = ind2depo(order, node)
    nbtn = length(pos_nat)
    dmax = maximum(depths)

    tmp = zeros(1, Int(2^dmax))
    beg = 1
    for k in 1:nbtn
        d = depths[k]
        len = 2^(dmax - d)
        tmp[Int(beg):Int(beg + len - 1)] .= d
        beg = beg + len
    end
    depths = tmp

    pos = Float64[]
    push!(pos, 0)
    for d in 1:dmax
        push!(pos, ((2^d - 1) .- pos)...)
    end
    spos = sort(pos) # uso una variabile tmp spos
    pos = indexin(spos, pos) # per ritornare gli indici nella funzione sort

    depths = depths[pos]
    pos = pos .+ 2^dmax .- 2

    for d in dmax-1:-1:1
        tmp = findall(depths .== d) #corrispettivo di find in matlab
        if !isempty(tmp)
            # da verificare!!!
            dd = dmax - d
            pow = 2^dd
            beg = tmp[1:Int(pow):end]
            tmp[1:Int(pow):end] = NaN
            pos[Int(beg)] = floor((pos[Int(beg)] + 1 - pow) / pow)
            pos[Int(tmp)] = NaN
        end
    end

    pos = setdiff(pos, NaN) # elimina i valori NaN e ricostruisce il vettore

    snode = sort(node)
    tmp = indexin(snode, node)
    spos = sort(pos)
    pos = indexin(spos, pos)
    spos = sort(pos)
    pos = indexin(spos, pos)

    ord = tmp[pos]
end

function otnodes(
    t::waveletPacketTree
)

    termnodes = []
    order = t.wavInfo.order
    node = t.tnodes.tn
    I = frqord(node)

    tn_Seq = node[I]
    sI = sort(I)
    J = indexin(sI, I)

    # if strcmpi(termnodes,"dp")
    #     [D,P] = ind2depo(ord,tn_Pal);
    #     tn_Pal = [D,P];
    #     [D,P] = ind2depo(ord,tn_Seq);
    #     tn_Seq = [D,P];    
    # end

    return node, tn_Seq, I, J
end

function wkeep1(
    z3::Union{AbstractVector{T},AbstractArray{T}}, # x
    len::Int64,
) where {T<:AbstractFloat}

    y = z3
    sx = length(z3)
    # ok = (len >= 0) && (len < sx) # sanity check
    side = 0
    d = (sx - len) / 2

    frst = Int(1 + floor(d))
    last = Int(sx - ceil(d))

    return y[frst:last]
end

function wpspectrum(
    x::Union{AbstractVector{T},AbstractArray{T}},
    sr::Int64,
    depth::Int64,
    wname::String;
    entType=:shannon::Symbol,
    entPar=0.0::Real
) where {T<:AbstractFloat}

    t = wpdec(x, depth, wname, entType, entPar)

    order = t.wavInfo.order
    dmax = t.wavInfo.depth
    tn = t.tnodes.tn
    sizes = t.tN.size
    nbtn = length(tn)
    cfs = t.tN.data

    depths, posis = ind2depo(order, tn)
    node, tn_Seq, I, J = otnodes(t)

    cfs = abs.(cfs)
    ord = I

    maxsizes = []
    for i in eachindex(sizes)
        push!(maxsizes, max(sizes[i]...)) #per riuscire a estrapolare il max da una tupla
    end
    sizes = maxsizes

    deb = ones(1, nbtn + 1)
    fin = ones(1, nbtn) # in matlab mette +1 ma secondo me è  un errore
    for k in 1:nbtn
        fin[k] = deb[k] + sizes[k] - 1
        deb[k+1] = fin[k] + 1
    end

    nbrows = (2 .^ (dmax .- depths))
    NBrowtot = sum(nbrows)
    datasize, channels = t.wavInfo.dataSize
    len = datasize
    spec = zeros(Int(NBrowtot), len)
    ypos = zeros(nbtn, 1)

    if nbtn > 1
        for k in 1:nbtn
            ypos[ord[k]] = sum(nbrows[ord[1:k-1]])
        end
    end

    ypos = NBrowtot .+ 1 .- ypos .- nbrows

    for k in 1:nbtn
        d = depths[k]
        z = cfs[Int.(deb[Int(k)]:fin[Int(k)])]
        # transpose(z) # inverte righe e colonne!

        z2 = fill(z, Int(2^d))
        r = size(z2)[1]
        c = size(z2[1])[1]
        z2 = reshape(reduce(vcat, z2), (c, r)) # trasforma da vettore di vettori, a matrice
        z3 = reshape(transpose(z2), 1, :)

        vals = wkeep1(z3, len)

        r1 = Int(ypos[k])
        r2 = Int(ypos[k] + nbrows[k] - 1)

        v1 = fill(vals, Int(nbrows[k]))
        r = size(v1)[1]
        c = size(v1[1])[1]
        spec[r1:r2, :] = reshape(reduce(vcat, v1), (c, r)) # trasforma da vettore di vettori, a matrice

        # if LowFreqOff && k == 1
        #     spec(r1:r2, :) = NaN
        # end
    end
    freqs=[]
    for i in 1:size(spec, 1)
        push!(freqs, i * sr / (2 * size(spec, 1)))
    end
    times = []
    for i in 0:1:datasize-1
        push!(times, i / sr)
    end
    times = (0:1:datasize-1) / sr

    return spec, freqs, times
end

# # debug
# using PyCall
# librosa = pyimport_conda("librosa")
# sr_src = 8000
# x, sr = librosa.load("/home/riccardopasini/Documents/Aclai/Julia_additional_files/test.wav", sr=sr_src, mono=true)

# depth = 5
# wname = "db2"
# spec, freqs, times = wpspectrum(x, sr, depth, wname)