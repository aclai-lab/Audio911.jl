# ---------------------------------------------------------------------------- #
#                     windows - adapted from DSP package                       #
# ---------------------------------------------------------------------------- #
function get_window(winfunc::Function, nwin::Int; padding::Int = 0, zerophase::Bool = false, periodic::Bool = false)
    @assert nwin > 0 "nwin must be positive"
    @assert padding ≥ 0 "padding must be nonnegative"    

    if nwin == 1
        # if nwin is set to 1, no windowing will be applied. for future applications with wavelets.
        return 1.0
    elseif zerophase
        return vcat([winfunc.(range(0.0, stop=(nwin÷2)/nwin, length=nwin÷2+1))], [winfunc.(range(-(nwin÷2)/nwin, stop=-1/nwin, length=nwin÷2))])
    elseif periodic
        return vcat(winfunc.(range(-0.5, stop=0.5, length=nwin+1))[1:end-1])
    else
        return vcat(winfunc.(range(-0.5, stop=0.5, length=nwin)))
    end
end

function hann(nwin::Int; kwargs...)
    get_window(nwin; kwargs...) do x
        0.5 * (1 + cospi(2x))
    end
end

function hamming(nwin::Int; kwargs...)
    get_window(nwin; kwargs...) do x
        0.54 - 0.46 * cospi(2x)
    end
end

function rect(nwin::Int; kwargs...)
    get_window(nwin; kwargs...) do _
        1.0
    end
end

### fix from here

# function tukey(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)
#     # check that α is reasonable
#     !(0 <= α <= 1) && error("α must be in the range 0 <= α <= 1.")

#     # if α is less than machine precision, call it zero and return the
#     # rectangular window for this length.  if we don't short circuit this
#     # here, it will blow up below.
#     abs(α) <= eps() && return rect(n; padding=padding, zerophase=zerophase)

#     makewindow(n, padding, zerophase) do x
#         if x <= -(1-α)/2
#             0.5*(1 + cos(2pi/α*(x+(1-α)/2)))
#         elseif x <= (1-α)/2
#             1.0
#         else
#             0.5*(1 + cos(2pi/α*(x-(1-α)/2)))
#         end
#     end
# end

# function cosine(n::Integer; padding::Integer=0, zerophase::Bool=false)
#     makewindow(n, padding, zerophase) do x
#         cos(pi*x)
#     end
# end

# function lanczos(n::Integer; padding::Integer=0, zerophase::Bool=false)
#     makewindow(n, padding, zerophase) do x
#         sinc(2x)
#     end
# end

# function triang(n::Integer; padding::Integer=0, zerophase::Bool=false)
#     # for the purpose of calculating the slope of the window, consider `n` to be
#     # 1 larger to compensate for the fact that `zerophase` gives a periodic
#     # window
#     m = zerophase ? n+1 : n
#     scale = iseven(m) ? 2(m-1)/m : 2(m-1)/(m+1)
#     makewindow(n, padding, zerophase) do x
#         1 - scale*abs(x)
#     end
# end

# function bartlett(n::Integer; padding::Integer=0, zerophase::Bool=false)
#     makewindow(n, padding, zerophase) do x
#         1 - abs(2x)
#     end
# end

# function gaussian(n::Integer, σ::Real; padding::Integer=0, zerophase::Bool=false)
#     σ > 0.0 || error("σ must be positive")
#     makewindow(n, padding, zerophase) do x
#         exp(-0.5*(x/σ)^2)
#     end
# end

# function bartlett_hann(n::Integer; padding::Integer=0, zerophase::Bool=false)
#     a0, a1, a2 = 0.62, 0.48, 0.38
#     makewindow(n, padding, zerophase) do x
#         a0 - a1*abs(x) + a2*cos(2pi*x)
#     end
# end

# function kaiser(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)
#     pf = 1.0/besseli(0,pi*α)
#     makewindow(n, padding, zerophase) do x
#         pf*besseli(0, pi*α*(sqrt(1 - (2x)^2)))
#     end
# end

# ---------------------------------------------------------------------------- #
#                                   windowing                                  #
# ---------------------------------------------------------------------------- #
function buffer(x::AbstractVector{T}, nwin::Int, noverlap::Int) where T<:AbstractFloat
    nhop = nwin - noverlap
    nhops = div(length(x) - nwin, nhop) + 1
    # hcat(map(i -> x[(i - 1) * nhop + 1 : (i - 1) * nhop + nwin], 1:nhops)...)

    y = Matrix{T}(undef, nwin, nhops)
    
    @inbounds @simd for i in 1:nhops
        y[:, i] = @view x[(i - 1) * nhop + 1 : (i - 1) * nhop + nwin]
    end

    return y
end

function _get_frames(x::AbstractVector{Float64}, window::Tuple{Symbol, Symbol}, nwin::Int64, noverlap::Int64)
    @assert window[2] in [:periodic, :symmetric] "window can be only :symmetric or :periodic"
	buffer(x, nwin, noverlap), getfield(@__MODULE__, window[1])(nwin; periodic=window[2] == :periodic)
end