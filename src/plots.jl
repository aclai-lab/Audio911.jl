# ---------------------------------------------------------------------------- #
#                              AbstractSpectrogram                             #
# ---------------------------------------------------------------------------- #
function Audio911.plot(s::AbstractSpectrogram; 
    db::Bool=true,
    freq_scale::Symbol=:log10,
    colormap::Symbol=:viridis,
    clims::Union{Nothing,Tuple}=nothing,
    kwargs...)
    
    spec = get_data(s)
    T    = eltype(spec)
    freq = collect(get_freq(s))
    sr   = get_sr(s)
    winsize = get_winsize(s)
    overlap = get_overlap(s)

    first_bin = findfirst(f -> f > 0, freq)
    spec = spec[first_bin:end, :]
    freq = freq[first_bin:end]
    
    hop_size = winsize - overlap
    nframes = size(spec, 2)
    time = (0:nframes-1) .* (hop_size / sr)
    
    plot_data = db ? 10 .* log10.(spec .+ eps(T)) : spec
    
    if clims === nothing
        clims = db ? (-80, maximum(plot_data)) : (0, maximum(plot_data))
    end
    
    p = Plots.heatmap(
        time, freq, plot_data;
        xlabel="Time (s)",
        ylabel="Frequency (Hz)",
        title="STFT Spectrogram",
        colorbar_title=db ? "Power (dB)" : "Power",
        color=colormap,
        clims=clims,
        yscale=freq_scale,
        kwargs...
    )
    
    # Format y-axis ticks to show actual Hz values
    if freq_scale in (:ln, :log10, :log2)
        fmin, fmax = extrema(freq)
        
        # Generate nice tick positions
        if freq_scale == :log10
            log_min = floor(Int, log10(fmin))
            log_max = ceil(Int, log10(fmax))
            tick_positions = [10.0^i for i in log_min:log_max]
            
            # Create labels with k for thousands
            tick_labels = map(tick_positions) do f
                if f >= 1000
                    string(Int(f ÷ 1000), "k")
                else
                    string(Int(f))
                end
            end
            
            Plots.yticks!(p, (tick_positions, tick_labels))
        end
    end
    
    return p
end