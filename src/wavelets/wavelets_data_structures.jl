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

struct dwtAttr{T<:Real}
    extMode::Symbol
    shift1D::Union{Int64,Float64}
    shift2D::AbstractVector{T}
end