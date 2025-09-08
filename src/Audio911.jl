module Audio911
using  Reexport

# using AudioReader: load
@reexport using AudioReader: @format_str, File, load

# using FFTW, DSP
# using LinearAlgebra
# using Parameters
# using SpecialFunctions
# using StatsBase
# using Statistics, Roots
# using NaNStatistics

# include("signalDataStructure.jl")
# include("audioFeaturesExtractor.jl")
# # windowing
# include("windowing/cswindows.jl")
# include("windowing/windows.jl")
# include("windowing/windowing.jl")
# # fft
# include("fft/conv.jl")
# include("fft/f0.jl")
# include("fft/fft.jl")
# include("fft/lin.jl")
# include("fft/mel.jl")
# include("fft/spectral.jl")
# # utils
# include("utils/histogram.jl")
# include("utils/speech_detector.jl")
# include("utils/in_out.jl")
# include("utils/trimaudio.jl")
# # wavelets
# include("wavelet/cwt.jl")
# # constant-q transform
# include("cqt/cqt.jl")

# # structures
# export AudioSetup, AudioData, AudioObj
# # audio features
# export audio_obj, get_features

# # utility functions
# export speech_detector
# export load_audio, save_audio, trim_audio, normalize_audio
# # wavelets
# export cwt, cwt_windowing
# get_fft!
# # TODO patch
# extractfeatures = 
# export extractfeatures

# # modular
# export get_frames, get_frames!, get_stft, get_stft!

# export get_frames2, get_stft2

end # module Audio911
