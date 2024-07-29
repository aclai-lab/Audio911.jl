# wavelet packet structures
# https://it.mathworks.com/help/wavelet/ref/wptree.html

# wavelet packet object structure functions
function getTerminalNodes( # get terminal nodes
    order::Int,
    depth::Int,
    spsch::AbstractVector{T}
) where {T<:Real}

    nbtn = order^depth
    if (order == 0)
        tn = 0
    elseif (order == 1)
        tn = depth
    else
        tn = collect((nbtn-1)/(order-1):(order*nbtn-order)/(order-1)) # genera un vettore da..a
    end

    if (order < 2) || (depth < 2)
        return # valutarne il ritorno
    end

    asc = zeros(nbtn, depth)
    asc[:, 1] .= tn
    for j = 1:depth-1
        asc[:, j+1] .= floor.((asc[:, j] .- 1) / order)
    end
    tab = Int.(asc .- order * floor.((asc .- 1) / order))
    for i in eachindex(tab)
        tab[i] = spsch[tab[i]]
    end

    icol = ones(nbtn, 1)

    for ic = depth:-1:2
        ir = filter(x -> x < ic, icol) # praticamente cerca tutti gli elementi < ic e l'output è un vettore. x -> è una funzione
        K = zeros(0)
        for i in eachindex(ir)
            (tab[i, ic] == 0 ? append!(K, ir[tab[i, ic]]) : nothing) # K contiene sono valori se tab è 0
        end
        K = Int.(K)
        for i in eachindex(K)
            icol[i] = ic
            tn[i] = asc[i, ic]
        end
    end
    ttn = zeros(0)
    append!(ttn, tn[1])
    for i in 1:length(tn)-1
        if (tn[i] - tn[i+1] != 0)
            append!(ttn, tn[i+1])
        end
    end
    return ttn, length(ttn)
end # function getTerminalNodes

function ntree(
    order=2::Int,
    depth=0::Int,
    spsch=Int.(ones(order))
)
    tn, nbtn = getTerminalNodes(order, depth, spsch)

    return tn, nbtn
end

function qmf(
    x::AbstractVector{T},
    p=0
) where {T<:Real}
    # quadrature mirror filter
    # Y = QMF(X,P) changes the signs of the even index entries of the reversed vector filter coefficients X if P is even.
    # If P is odd the same holds for odd index entries. Y = QMF(X) is equivalent to Y = QMF(X,0).

    y = reverse(x)
    first = 2 - rem(p, 2) # rem = remainder after division
    y[first:2:end] = -y[first:2:end]

    return y
end # function qmf

function orthfilt(
    w::AbstractVector{T},
    p=0 # se P non viene passato, va impostato a 0
) where {T<:Real}

    w = w / sum(w)

    lo_R = sqrt(2) .* w
    hi_R = qmf(lo_R, p)
    lo_D = reverse(lo_R)
    hi_D = reverse(hi_R)

    return lo_D, hi_D, lo_R, hi_R
end # function orthfilt 

function wfilters(
    wname::String
)
    file = match(r"(?:[a-zA-Z]+)", wname).match # estrapola le lettere tramite regex da wname, .match riconverte da regex a string
    wcode = wname[length(file)+1:end] # estrapola la parte numerica
    i_fam = winfo[file] # recupera i dati dal dizionario generale

    F = i_fam.coeff[wcode]
    lo_D, hi_D, lo_R, hi_R = orthfilt(F)
end # function wFilters

function wunique(
    x::AbstractVector
) # remove duplicates, return sorted array

    if (isempty(x))
        return y = []
    end

    y = sort(x)
    d = diff(y)
    J = 1 .+ findall(x -> x > 0, d) # ritorna gli indici degli elementi >0
    K = [1; J]
    y = y[K]
end # function wunique

function ascendants(
    nodes::Vector{Float64},
    order::Int64,
    level::Int64,
    flagdp::Union{Symbol,Bool}
) # construction of ascendants table

    tab = zeros(length(nodes), level + 1)
    tab[:, 1] = sort(nodes)
    for j = 1:level
        tab[:, j+1] = floor.((tab[:, j] .- 1) / order)
    end

    idx = [ones.(1, level + 1); diff(tab, dims=1)]
    # diff trova di quanto variano i valori lungo una dimensione
    # in idx[... ;...] si sommano appendono 2 matrici
    tab = tab[idx.>0] # applica idx come maschera
    tab = wunique(tab) # tab = wunique(tab[:])
    tab = tab[tab.>=0]
    # if nargin==4 && flagdp
    #     [tab(:,1),tab(:,2)] = ind2depo(order,tab);
    # end 
end # function ascendants

function allnodes(
    t::waveletPacketTree
)
    order = t.wavInfo.order
    depth = t.wavInfo.depth
    allN = t.tnodes.tn

    if (depth == 0)
        return
    end
    flagdp = :false

    allN = ascendants(allN, order, depth, flagdp)
end # function all nodes

function depo2ind(
    order::Int64,
    nodes::AbstractVector
) # node depth-position to node index
    r = size(nodes, 1)
    c = size(nodes, 2)
    if (c == 0 || c == 1)
        return nodes
    elseif (c == 2)
        # da vedere
    else
        #error
    end
end # function depo2ind

function ismemberBuiltinTypes(
    a::AbstractArray{T},
    b::AbstractArray{T}
) where {T<:Real}
    locb = zeros(size(a))
    numelA = length(a)
    numelB = length(b)

    # if (numelA == 0 || numelB <= 1)
    #     if (numelA > 0 && numelB == 1)
    #         lia = (a == b)
    #         locb = double(lia)
    #     else
    #         lia = false(size(a))
    #     end
    #     return
    # end 

    scalarcut = 5
    if (numelA <= scalarcut)
        # lia = false(size(a))
        # for i=1:numelA
        #     found = a(i)==b(:);
        #     if any(found)
        #         lia(i) = true;
        #         locb(i) = find(found,1);
        #     end
        # end
    else
        sortedlist = issorted(b[:])
        if (!sortedlist)
            idx = sortperm(b[:])
            b = sort(b[:])
        end

        lia = map(in(b), a)
        locb = zeros(size(a))
        for i in eachindex(a)
            if (a[i] in (b))
                locb[i] = findfirst(x -> x == a[i], b)
            end
        end
    end

    return lia, locb
end # function ismemberBuiltinTypes

function ismember(
    a::AbstractArray{T},
    b::AbstractArray{T}
) where {T<:Real} # for arrays A and B returns an array of the same size as A containing true where the elements of A are in B and false otherwise

    byrow = :false
    doBuiltinTypes = :true

    # check that one of A and B is double if A and B are non-homogeneous. 
    # do a separate check if A is a heterogeneous object and only allow a B that is of the same root class

    if (!byrow)
        if (doBuiltinTypes)
            lia, locb = ismemberBuiltinTypes(a, b)
        else
            # lia,locb = ismemberClassTypes(a,b)
        end
    else
        # row case
    end

    return lia, locb
end

function istnode(
    t::waveletPacketTree,
    n::AbstractVector
) # determine indices of terminal nodes
    order = t.wavInfo.order
    tn = t.tnodes.tn
    n = depo2ind(order, n)
    lia, r = ismember(n, tn)

    return r
end # function istnode

function locnumcn(
    order::Int64,
    nodes::AbstractVector{T}
) where {T<:Real} # local number for a child node
    if (order == 0)
        return 1
    end
    return nodes .- order * floor.((nodes .- 1) / order)
end # function locnumcn

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

function dtree(
    x::Union{AbstractVector{T},AbstractArray{T}},
    order=2::Int,
    depth=0::Int,
    spflg=:notexpand::Symbol,
    spsch=Int.(ones(order))
) where {T<:Real}

    tn, nbtn = ntree(order, depth, spsch)

    tn = Int.(tn)

    if (spflg == :expand)
        # EXPAND Expand data tree.
        # expand.m
        # da implementare
    end

    d = Vector{allNodes}(undef, tn[end] + 1)
    d[1] = allNodes(0, (size(x, 1), size(x, 2)), NaN, NaN)
    for i in 1:tn[end]
        d[i+1] = allNodes(i, (0, 0), NaN, NaN)
    end

    return wTnodes(tn, nbtn, spsch), d
end

function wptree(
    x::Union{AbstractVector{T},AbstractArray{T}},
    wname::String,
    order=2::Int64,
    depth=0::Int64,
    entType=:shannon::Symbol,
    entPar=0.0::Real
) where {T<:Real}

    wtn, d = dtree(x, order, depth, :notexpand)

    lo_D, hi_D, lo_R, hi_R = wfilters(wname)

    winfo = wavInfo(
        x, (size(x, 1), size(x, 2)),
        order, depth,
        wname,
        lo_D, hi_D, lo_R, hi_R
    )

    went = entInfo(entType, entPar)

    n=terminalNodes([], [0 0; 0 0])

    t = waveletPacketTree(d, n, wtn, winfo, went)
    t = expand(t)
end # function wptree

# MAIN avelet packet decomposition main function
function wpdec(
    x::Union{AbstractVector{T},AbstractArray{T}},
    depth::Int64,
    wname::String,
    entType=:shannon::Symbol,
    entPar=0.0::Real
) where {T<:Real}
    ## TODO check args

    order = 2

    t = wptree(x, wname, order, depth, entType, entPar) # tree computation

end # function wpdec

