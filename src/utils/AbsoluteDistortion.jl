## main
using PyCall
librosa = pyimport_conda("librosa")
soundfile = pyimport_conda("soundfile")

sr_src = 48000
x, sr = librosa.load("/media/paso/Paso Stick/Trombone Live_01.wav", sr=sr_src, mono=true)

absX = abs.(x)

soundfile.write("test1.wav",absX,samplerate=sr,subtype="PCM_24")