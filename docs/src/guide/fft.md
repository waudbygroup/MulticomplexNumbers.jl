# [FFT Operations](@id fft-operations)

```@meta
CurrentModule = MulticomplexNumbers
```

This page explains how to perform Fourier transforms on multicomplex data, particularly for NMR spectroscopy applications.

---

## Setup

FFT support requires the FFTW package:

```julia
using MulticomplexNumbers
using FFTW  # This loads the FFT extension
```

---

## Basic FFT

The `fft!` function performs in-place FFT along a specific imaginary unit:

```julia
using MulticomplexNumbers
using FFTW

# Create bicomplex data
data = rand(Multicomplex{Float64,2,4}, 64)

# FFT along im1 (treating im1 as the complex unit)
fft!(data, 1)

# FFT along im2
fft!(data, 2)
```

---

## Unit-Dimension Association

**This is the most common pattern for NMR data processing.**

In multi-dimensional NMR, we typically associate each imaginary unit with a specific array dimension:

```julia
using MulticomplexNumbers
using FFTW

# Create 2D NMR data (64 × 32 bicomplex array)
# Dimension 1 (rows): direct dimension, associated with im1
# Dimension 2 (cols): indirect dimension, associated with im2
data = Array{Multicomplex{Float64, 2, 4}}(undef, 64, 32)

# ... fill with time-domain data ...

# Transform the direct dimension (im1) along dimension 1
fft!(data, 1, dims=1)

# Transform the indirect dimension (im2) along dimension 2
fft!(data, 2, dims=2)
```

This pattern extends naturally to higher dimensions:

```julia
using MulticomplexNumbers
using FFTW

# 3D NMR data (tricomplex)
data_3d = Array{Multicomplex{Float64, 3, 8}}(undef, 128, 64, 32)

# ... fill with time-domain data ...

# Transform each dimension with its associated imaginary unit
fft!(data_3d, 1, dims=1)  # Direct dimension (im1)
fft!(data_3d, 2, dims=2)  # First indirect (im2)
fft!(data_3d, 3, dims=3)  # Second indirect (im3)
```

---

## Custom FFT Dimensions

You can apply FFT along specific array dimensions:

```julia
using MulticomplexNumbers
using FFTW

# Create 3D array with bicomplex elements
data = Array{Multicomplex{Float64, 2, 4}}(undef, 32, 32, 32)

# ... fill with data ...

# FFT along im1, but only along dimension 2 of the array
fft!(data, 1, dims=2)

# FFT along im2, along both dimensions 1 and 3
fft!(data, 2, dims=(1, 3))

# Most common: one unit per dimension
fft!(data, 1, dims=1)  # im1 transforms along array dim 1
fft!(data, 2, dims=2)  # im2 transforms along array dim 2
# Dimension 3 could be e.g., different experiments
```

---

## Inverse FFT

Use `ifft!` to perform the inverse Fourier transform. This is essential for NMR reprocessing workflows (e.g. solvent suppression, linear prediction):

```julia
using MulticomplexNumbers
using FFTW

data = Array{Multicomplex{Float64, 2, 4}}(undef, 64, 32)
# ... fill with time-domain data ...

# Forward transform
fft!(data, 1, 1)
fft!(data, 2, 2)

# ... apply corrections in frequency domain ...

# Inverse transform back to time domain
ifft!(data, 1, 1)
ifft!(data, 2, 2)
```

An unnormalized inverse FFT is also available via `bfft!`, where `bfft!(fft!(copy(A), n), n) ≈ A .* length(A)`.

---

## Allocating FFT

If you need to preserve the original data, use the allocating variants `fft` and `ifft` (without the `!`):

```julia
spectrum = fft(fid, 1)    # fid is unchanged
fid_back = ifft(spectrum, 1)  # spectrum is unchanged
```

---

## Multicomplex Unit Dispatch

You can pass a multicomplex imaginary constant directly instead of an integer:

```julia
fft!(data, im1, 1)  # equivalent to fft!(data, 1, 1)
fft!(data, im2, 2)  # equivalent to fft!(data, 2, 2)
```

---

## Programmatic Imaginary Units

For loops over dimensions, use `imN(n)` to construct the n-th imaginary unit:

```julia
# Process all dimensions of an N-dimensional experiment
for n in 1:N
    fft!(data, n, n)
    data .*= exp(imN(n) * phase[n])  # phase correct each dimension
end
```

---

## Frequency Shifting

`fftshift` and `ifftshift` from FFTW/AbstractFFTs work directly on multicomplex arrays:

```julia
using MulticomplexNumbers
using FFTW

fid = rand(Multicomplex{Float64,2,4}, 64)
fft!(fid, 1)

# Center the zero-frequency component
spectrum = fftshift(fid)
```

---

## See Also

- **[NMR Spectroscopy Applications](@ref nmr-applications)**: Complete NMR workflow examples
- **[Accessing Components](@ref components)**: Extract parts and components
- **[API Reference](@ref)**: Complete function documentation
