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
data = [Multicomplex(rand(4)...) for _ in 1:64]

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

## Complete Example: 2D NMR Processing

```julia
using MulticomplexNumbers
using FFTW

# Create synthetic 2D NMR signal
function create_2d_nmr_signal(n1, n2, freq1, freq2)
    signal = Array{Multicomplex{Float64, 2, 4}}(undef, n1, n2)

    for j in 1:n2, i in 1:n1
        # Time points
        t1 = (i-1) / n1
        t2 = (j-1) / n2

        # Oscillating signal with decay
        # Direct dimension (im1)
        s1_real = cos(2π * freq1 * t1) * exp(-t1 * 5.0)
        s1_imag = sin(2π * freq1 * t1) * exp(-t1 * 5.0)

        # Indirect dimension (im2)
        s2_real = cos(2π * freq2 * t2) * exp(-t2 * 3.0)
        s2_imag = sin(2π * freq2 * t2) * exp(-t2 * 3.0)

        # Combine as bicomplex
        signal[i, j] = Multicomplex(
            s1_real * s2_real,     # real-real
            s1_imag * s2_real,     # im1 coefficient
            s1_real * s2_imag,     # im2 coefficient
            s1_imag * s2_imag      # im1*im2 coefficient
        )
    end

    return signal
end

# Generate signal
fid = create_2d_nmr_signal(128, 64, 10.0, 15.0)

# Process: FFT both dimensions
fft!(fid, 1, dims=1)  # Transform direct dimension
fft!(fid, 2, dims=2)  # Transform indirect dimension

# The peak should appear near indices corresponding to the frequencies
# (11, 16) for frequencies 10 and 15 out of 128 and 64 points
```

---

## Phase Correction

Phase correction is performed by multiplying by complex exponentials. For multicomplex numbers, we use different imaginary units for different dimensions:

```julia
using MulticomplexNumbers
using FFTW

# After FFT, you may need to apply phase correction

# Zero-order phase correction for the direct dimension (im1)
function phase_correct_direct!(data, theta1)
    phase_factor = exp(im1 * theta1)
    data .*= phase_factor
    return data
end

# Zero-order phase correction for the indirect dimension (im2)
function phase_correct_indirect!(data, theta2)
    phase_factor = exp(im2 * theta2)
    data .*= phase_factor
    return data
end

# First-order phase correction (linear phase)
function phase_correct_first_order!(data, theta0, theta1_per_point, dim)
    if dim == 1
        n = size(data, 1)
        for i in 1:n
            phase = theta0 + theta1_per_point * (i - 1)
            data[i, :] .*= exp(im1 * phase)
        end
    elseif dim == 2
        n = size(data, 2)
        for j in 1:n
            phase = theta0 + theta1_per_point * (j - 1)
            data[:, j] .*= exp(im2 * phase)
        end
    end
    return data
end

# Example usage
fid = create_2d_nmr_signal(64, 32, 5.0, 10.0)
fft!(fid, 1, dims=1)
fft!(fid, 2, dims=2)

# Apply zero-order phase corrections
phase_correct_direct!(fid, π/6)      # 30° in direct dimension
phase_correct_indirect!(fid, π/4)    # 45° in indirect dimension

# Or first-order correction
# phase_correct_first_order!(fid, 0.0, π/180, 1)  # 1° per point in dim 1
```

---

## Extracting Complex Sub-Components

After processing, you often want to extract specific complex projections:

```julia
using MulticomplexNumbers
using FFTW

# Process 2D spectrum
fid = create_2d_nmr_signal(64, 32, 5.0, 10.0)
fft!(fid, 1, dims=1)
fft!(fid, 2, dims=2)

# Extract the real/im1 complex spectrum (traditional F1 dimension)
# This gives you the spectrum with im1 as the complex unit
spectrum_im1 = real.(fid)  # Returns Multicomplex{Float64, 1, 2} array

# Extract as actual Complex numbers for plotting/analysis
complex_im1 = ascomplex([spectrum_im1[i, j] for i in 1:64, j in 1:32], 1)

# Extract the real/im2 complex spectrum
# For each point, we want real + im2*i
function extract_im2_component(data)
    n1, n2 = size(data)
    result = Array{Complex{Float64}}(undef, n1, n2)
    for j in 1:n2, i in 1:n1
        z = data[i, j]
        result[i, j] = Complex(realest(z), component(z, 3))
    end
    return result
end

complex_im2 = extract_im2_component(fid)

# Or use ascomplex for the full array
# This creates a view that maps the specified unit to Complex{T}
complex_view = ascomplex(fid, 1)  # Map im1 → im (complex unit)
```

---

## Practical Tips for NMR

1. **Dimension-Unit Association**: Always associate `im1` with the direct (acquisition) dimension, `im2` with the first indirect, `im3` with the second indirect, etc.

2. **Processing Order**: Process dimensions in the order that makes sense for your experiment, typically:
   ```julia
   fft!(data, 1, dims=1)  # Direct first
   fft!(data, 2, dims=2)  # Then indirect
   ```

3. **Phase Correction**: Apply phase corrections using the appropriate imaginary unit for each dimension:
   - Direct dimension: multiply by `exp(im1 * θ)`
   - First indirect: multiply by `exp(im2 * θ)`
   - Second indirect: multiply by `exp(im3 * θ)`

4. **Memory Efficiency**: FFT is performed in-place, modifying the original array. Make a copy if you need to preserve the original:
   ```julia
   fid_copy = copy(fid)
   fft!(fid_copy, 1, dims=1)
   ```

5. **Extracting Real Spectra**: For publication, you typically want the magnitude:
   ```julia
   magnitude_spectrum = abs.(realest.(fid))
   ```

---

## Advanced: Custom FFT Dimensions

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

## See Also

- **[NMR Spectroscopy Applications](@ref nmr-applications)**: Complete NMR workflow examples
- **[Accessing Components](@ref components)**: Extract parts and components
- **[API Reference](@ref)**: Complete function documentation
