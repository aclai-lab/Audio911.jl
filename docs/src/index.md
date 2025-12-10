```@meta
CurrentModule = Audio911
```

```@raw html
<div align="center">
    <img src="../img/logo.png" alt="Audio911" width="600">
</div>

<h2 align="center">Audio Feature Extraction in Julia
<p align="center">
  <a href="https://github.com/aclai-lab/Audio911.jl/actions">
    <img src="https://github.com/aclai-lab/Audio911.jl/workflows/CI/badge.svg"
         alt="Build Status">
  </a>
  <a href="https://aclai-lab.github.io/Audio911.jl/dev/">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg"
         alt="dev documentation">
  </a>
  <a href="https://aclai-lab.github.io/Audio911.jl/stable/">
    <img src="https://img.shields.io/badge/docs-stable-blue.svg"
         alt="stable documentation">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yelllow"
       alt="bibtex">
  </a>
  <a href="https://codecov.io/gh/aclai-lab/Audio911.jl">
    <img src="https://codecov.io/gh/aclai-lab/Audio911.jl/branch/main/graph/badge.svg"
       alt="CodeCov">
  </a>
</p>
</h2>
```

**Audio911.jl** is your Swiss Army knife for extracting audio features in a simple and fast way.

Inspired by MATLAB's audio feature extraction toolkit, Audio911.jl guarantees the same results while being designed to be modular, allowing you to connect various extraction algorithms as you prefer. It currently provides STFT, linear spectrograms, mel/bark/ERB spectrograms, and MFCC coefficients, with new algorithms being added constantly—so stay tuned!

## Features

### Time-Frequency Representations
- **STFT**: Short-Time Fourier Transform with customizable windows
- **Linear Spectrogram**: `LinSpec`
- **Mel Spectrogram**: `MelSpec` (HTK and Slaney styles)
- **Bark Spectrogram**: Auditory Bark scale
- **ERB Spectrogram**: Equivalent Rectangular Bandwidth scale

### Coefficients
- **MFCC**: Mel-Frequency Cepstral Coefficients with delta/delta-delta
- **Customizable Rectification**: Log or cubic root
- **Energy Options**: Standard or MFCC-based log energy

### Extras
- **Modular Design**: Compose your own audio processing pipelines by chaining algorithms
- **Multi-Format Audio**: Load WAV, MP3, FLAC and OGG files via [AudioReader.jl](https://github.com/PasoStudio73/AudioReader.jl)
- **On-the-Fly Resampling**: Built-in sample rate conversion

## Installation

```julia
using Pkg
Pkg.add("Audio911")
```

## Quick Start

```julia
using Audio911

# Load an audio file (automatically resampled to 16kHz)
audio = load("speech.wav"; sr=16000, mono=true, norm=true)

# Compute STFT
stft = Stft(audiofile; win=movingwindow(winsize=512, winstep=256), type=hamming, periodic=true, spectrum=power)

# Generate Mel spectrogram
mel_spec = MelSpec(stft; win_norm=true, nbands=32, norm=bandwidth, domain=:linear, scale=htk)

# Extract MFCC coefficients
mfcc = Mfcc(mel_spec; ncoeffs=30, rect=cubic_root)

# Access the results
mfcc_coeffs = get_data(mfcc)  # 13×N matrix of coefficients
```

## Learn More

For a comprehensive understanding of Audio911.jl's capabilities, check out our **[tutorials in the documentation](https://aclai-lab/Audio911.jl/stable/tutorials/)**. The tutorials cover:

- **Step-by-step guide** for building complete audio processing pipelines
- **Real-world examples** for speech recognition, music analysis, and environmental sound classification
- **Best practices** for parameter selection and optimization
- **Advanced techniques** including custom filterbank designs and feature engineering

Each tutorial includes reproducible code examples with real audio files, so you can follow along and adapt the techniques to your own projects.

## About

Audio911.jl is developed by the [ACLAI Lab](https://aclai.unife.it/en/) @ University of Ferrara.

## License

MIT License