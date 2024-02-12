# reference: ETSI ES 201 108 V1.1.2 (2000-04)
# https://www.3gpp.org/ftp/tsg_sa/TSG_SA/TSGS_13/Docs/PDF/SP-010566.pdf

# use FFTW package

# function get_fft_length( # switch case da ricordare!!!
#     sr::Int64,
# )
#     sr <= 4000 && return 128
#     sr <= 8000 && return 256
#     sr <= 16000 && return 512
#     return 1024
# end # get_fft_length

function get_onesided_fft_range(fft_length::Int64)
    if mod(fft_length, 2) == 0
        return collect(1:Int(fft_length / 2 + 1))   # EVEN
    else
        return collect(1:Int((fft_length + 1) / 2))  # ODD
    end
end # get_onesided_fft_range

function takeFFT(
    data::signal_data,
    setup::signal_setup
)
    # TODO validate FFT length, validate overlap length (audioFeatureExtractor.m line 1324)

    setup.fft_length = setup.window_length # definisce la fft pari alla finestra
    hop_length = setup.window_length - setup.overlap_length
    data.fft_window, unused = gencoswin(setup.window_type[1], setup.window_length, setup.window_type[2])

    # split in windows
    y = buffer(data.x, setup.window_length, hop_length)

    # apply window and take fft
    Z = fft(y .* data.fft_window, (1,))

    # take one side
    logical_ossb = falses(setup.fft_length)
    logical_ossb[get_onesided_fft_range(setup.fft_length)] .= true
    Z = Z[logical_ossb, :]

    setup.spectrum_type == :power ? data.fft = real(Z .* conj(Z)) : data.fft = abs.(Z)

    # log energy
    # reference: ETSI ES 201 108 V1.1.2 (2000-04)
    # https://www.3gpp.org/ftp/tsg_sa/TSG_SA/TSGS_13/Docs/PDF/SP-010566.pdf
    if setup.log_energy_source == :standard
        log_energy = sum(eachrow(y .^ 2))

        if setup.normalization_type == :standard
            log_energy[log_energy.==0] .= floatmin(Float64)
        elseif setup.normalization_type == :dithered
            log_energy[log_energy.<1e-8] .= 1e-8
        end
        
        data.log_energy = log.(log_energy)
    end

end # takeFFT(data, setup)

function takeFFT(
    x::AbstractArray{T},
    sr::Int64;
    fft_length::Int64=256,
    window_type::Vector{Symbol}=[:hann, :periodic],
    window_length::Int64=Int(round(0.03 * sr)),
    overlap_length::Int64=Int(round(0.02 * sr)),
    window_norm::Bool=true,
    frequency_range::Vector{Int64}=[0, Int(floor(sr /2))],
    spectrum_type::Bool=:power
) where {T<:AbstractFloat}
    # setup and data structures definition    
    setup = signal_Setup(
        sr=sr,
        fft_length=fft_length,
        window_type=window_type,
        window_length=window_length,
        overlap_length=overlap_length,
        window_norm=window_norm,
        frequency_range=frequency_range,
        spectrum_type=spectrum_type
    )

    data = signal_data(
        x=Float64.(x)
    )

    takeFFT(data, setup)

    return data.fft, setup.fft_frequencies
end # takeFFT(kwarg...)