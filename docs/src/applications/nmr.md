# [NMR Spectroscopy](@id nmr-applications)

```@meta
CurrentModule = MulticomplexNumbers
```

This tutorial demonstrates how to use multicomplex numbers for multi-dimensional NMR spectroscopy data processing.

---

## Background: Why Multicomplex Numbers for NMR?

In NMR spectroscopy:
- Each time-domain dimension uses **quadrature detection**, recording both real and imaginary components
- A 2D NMR experiment produces data that is complex in both the direct and indirect dimensions
- Traditional approaches store this as complex arrays with special conventions for higher dimensions

**Multicomplex numbers provide a natural, unified representation:**
- Order 2 (bicomplex) for 2D NMR: `z = a + b*i₁ + c*i₂ + d*i₁i₂`
- Order 3 (tricomplex) for 3D NMR: 8 components representing all combinations
- Order 4 for 4D NMR, etc.

Where:
- `i₁` represents the imaginary component of the direct (acquisition) dimension
- `i₂` represents the imaginary component of the first indirect dimension
- `i₃` represents the imaginary component of the second indirect dimension
- And so on...

---

## Tutorial: Processing 2D NMR Data

### Step 1: Setup

```julia
using MulticomplexNumbers
using FFTW  # Loads the FFT extension
```

### Step 2: Creating Synthetic NMR Data

Let's create a synthetic 2D NMR signal to demonstrate the workflow:

```julia
"""
Create a synthetic 2D NMR time-domain signal.

# Arguments
- `n1::Int`: Number of points in direct dimension
- `n2::Int`: Number of points in indirect dimension
- `freq1::Float64`: Frequency in direct dimension (Hz)
- `freq2::Float64`: Frequency in indirect dimension (Hz)
- `decay1::Float64`: Decay rate in direct dimension (1/s)
- `decay2::Float64`: Decay rate in indirect dimension (1/s)
"""
function create_2d_fid(n1, n2, freq1, freq2, decay1=5.0, decay2=3.0)
    # Allocate bicomplex array
    fid = Array{Multicomplex{Float64, 2, 4}}(undef, n1, n2)

    # Generate time-domain signal
    for j in 1:n2
        t2 = (j-1) / n2  # Indirect dimension time

        # Indirect dimension signal
        s2_r = cos(2π * freq2 * t2) * exp(-t2 * decay2)
        s2_i = sin(2π * freq2 * t2) * exp(-t2 * decay2)

        for i in 1:n1
            t1 = (i-1) / n1  # Direct dimension time

            # Direct dimension signal
            s1_r = cos(2π * freq1 * t1) * exp(-t1 * decay1)
            s1_i = sin(2π * freq1 * t1) * exp(-t1 * decay1)

            # Combine as bicomplex: (s1_r + s1_i*im1) * (s2_r + s2_i*im2)
            fid[i, j] = Multicomplex(
                s1_r * s2_r,     # Real-real component
                s1_i * s2_r,     # im1 coefficient
                s1_r * s2_i,     # im2 coefficient
                s1_i * s2_i      # im1*im2 coefficient
            )
        end
    end

    return fid
end

# Create test data: 128 × 64 points, peaks at 10 Hz (F1) and 15 Hz (F2)
fid = create_2d_fid(128, 64, 10.0, 15.0)
```

### Step 3: Fourier Transform

**Key principle:** Associate each imaginary unit with its corresponding array dimension.

```julia
# Transform direct dimension (im1 along dimension 1)
fft!(fid, 1, dims=1)

# Transform indirect dimension (im2 along dimension 2)
fft!(fid, 2, dims=2)

# Now fid contains the 2D frequency-domain spectrum
```

This pattern is clean and intuitive: `im1` ↔ dimension 1, `im2` ↔ dimension 2.

### Step 4: Phase Correction

Phase errors are common in NMR. We correct them by multiplying by phase factors.

#### Zero-Order Phase Correction

```julia
"""Apply zero-order phase correction to direct dimension."""
function phase_correct_0th_direct!(spectrum, theta1)
    spectrum .*= exp(im1 * theta1)
    return spectrum
end

"""Apply zero-order phase correction to indirect dimension."""
function phase_correct_0th_indirect!(spectrum, theta2)
    spectrum .*= exp(im2 * theta2)
    return spectrum
end

# Example: Apply 30° correction to direct, 45° to indirect
phase_correct_0th_direct!(fid, π/6)      # 30° = π/6 radians
phase_correct_0th_indirect!(fid, π/4)    # 45° = π/4 radians
```

#### First-Order Phase Correction

First-order phase correction applies a linearly varying phase:

```julia
"""
Apply first-order phase correction.

# Arguments
- `spectrum`: The spectrum array
- `theta0`: Zero-order phase (radians)
- `theta1`: First-order phase (radians per point)
- `dim`: Dimension to correct (1 or 2)
"""
function phase_correct_1st!(spectrum, theta0, theta1, dim)
    if dim == 1
        n = size(spectrum, 1)
        for i in 1:n
            phase = theta0 + theta1 * (i - 1) / n
            spectrum[i, :] .*= exp(im1 * phase)
        end
    elseif dim == 2
        n = size(spectrum, 2)
        for j in 1:n
            phase = theta0 + theta1 * (j - 1) / n
            spectrum[:, j] .*= exp(im2 * phase)
        end
    end
    return spectrum
end

# Apply first-order correction: 0° to 180° across direct dimension
phase_correct_1st!(fid, 0.0, π, 1)
```

### Step 5: Extracting Complex Sub-Components

After processing, you often need to extract specific complex projections for visualization or further analysis.

#### Extracting Real/Im1 Component

This gives you the spectrum as if only `im1` were the complex unit (traditional F1 view):

```julia
# Method 1: Using real() to get order-1 multicomplex
spectrum_im1 = real.(fid)  # Returns array of Multicomplex{Float64, 1, 2}

# Convert to standard Complex for plotting
complex_f1 = [Complex(realest(z), imag(z)) for z in spectrum_im1]

# Method 2: Using ascomplex (creates a view)
# Note: ascomplex expects contiguous memory, so we may need to reshape
spectrum_1d = vec(spectrum_im1)
complex_view = ascomplex(spectrum_1d, 1)
```

#### Extracting Real/Im2 Component

This gives you the spectrum with `im2` as the complex unit:

```julia
"""Extract the real + im2*i component from bicomplex data."""
function extract_real_im2(data::Array{<:Multicomplex{T,2,4}}) where T
    n1, n2 = size(data)
    result = Array{Complex{T}}(undef, n1, n2)

    for j in 1:n2, i in 1:n1
        z = data[i, j]
        # Component 1 is real, component 3 is im2 coefficient
        result[i, j] = Complex(component(z, 1), component(z, 3))
    end

    return result
end

complex_f2 = extract_real_im2(fid)
```

#### Extracting Magnitude Spectrum

For publication, you typically want the magnitude (absolute value):

```julia
# Most real component magnitude
magnitude = abs.(realest.(fid))

# Or full multicomplex magnitude
magnitude_full = abs.(fid)
```

---

## Complete 2D NMR Processing Example

Here's a complete workflow:

```julia
using MulticomplexNumbers
using FFTW

# 1. Create or load data
fid = create_2d_fid(256, 128, 12.5, 8.3)

# 2. Apply apodization (window function) if desired
function apply_exponential_window!(data, lb1, lb2)
    n1, n2 = size(data)

    for j in 1:n2, i in 1:n1
        t1 = (i-1) / n1
        t2 = (j-1) / n2
        window = exp(-π * lb1 * t1) * exp(-π * lb2 * t2)
        data[i, j] *= window
    end

    return data
end

apply_exponential_window!(fid, 2.0, 1.0)  # 2 Hz in F1, 1 Hz in F2

# 3. Zero-fill (optional, for higher digital resolution)
function zero_fill(data, new_n1, new_n2)
    n1, n2 = size(data)
    T = eltype(data)
    result = zeros(T, new_n1, new_n2)
    result[1:n1, 1:n2] .= data
    return result
end

fid = zero_fill(fid, 512, 256)

# 4. FFT
fft!(fid, 1, dims=1)  # Direct dimension
fft!(fid, 2, dims=2)  # Indirect dimension

# 5. Phase correction
phase_correct_0th_direct!(fid, π/8)      # 22.5°
phase_correct_0th_indirect!(fid, π/6)    # 30°
phase_correct_1st!(fid, 0.0, π/4, 1)     # 45° across F1

# 6. Extract magnitude for plotting
magnitude = abs.(realest.(fid))

# 7. Find peak maximum
max_val, max_idx = findmax(magnitude)
println("Peak at index $(max_idx.I) with intensity $(max_val)")

# 8. Extract complex projections if needed
f1_projection = extract_real_im2(fid)    # F1 dimension
f2_projection = [Complex(realest(z), imag(z)) for z in real.(fid)]  # F2 dimension
```

---

## 3D NMR Processing

For 3D experiments, use tricomplex numbers (order 3):

```julia
using MulticomplexNumbers
using FFTW

# Create 3D data
fid_3d = Array{Multicomplex{Float64, 3, 8}}(undef, 128, 64, 32)

# ... fill with experimental data ...

# Transform each dimension with its associated unit
fft!(fid_3d, 1, dims=1)  # Direct dimension (im1)
fft!(fid_3d, 2, dims=2)  # First indirect (im2)
fft!(fid_3d, 3, dims=3)  # Second indirect (im3)

# Phase correction for each dimension
fid_3d .*= exp(im1 * 0.1)   # Direct
fid_3d .*= exp(im2 * 0.2)   # First indirect
fid_3d .*= exp(im3 * 0.15)  # Second indirect

# Extract magnitude
magnitude_3d = abs.(realest.(fid_3d))
```

---

## Advanced: Selective Dimension Processing

Sometimes you want to process only specific dimensions:

```julia
using MulticomplexNumbers
using FFTW

# 2D array of bicomplex numbers
data = create_2d_fid(128, 128, 10.0, 15.0)

# Process only the direct dimension
fft!(data, 1, dims=1)
# Now: frequency in F1, time in F2

# Apply phase correction only to the transformed dimension
data .*= exp(im1 * π/4)

# Later, transform F2
fft!(data, 2, dims=2)
data .*= exp(im2 * π/6)
```

---

## Tips for Real NMR Data

1. **Data Loading**: Convert your NMR data to bicomplex/tricomplex format:
   ```julia
   function nmr_to_bicomplex(real_part, imag_part_1, imag_part_2, cross_part)
       n1, n2 = size(real_part)
       result = Array{Multicomplex{Float64, 2, 4}}(undef, n1, n2)

       for j in 1:n2, i in 1:n1
           result[i, j] = Multicomplex(
               real_part[i, j],
               imag_part_1[i, j],
               imag_part_2[i, j],
               cross_part[i, j]
           )
       end

       return result
   end
   ```

2. **States-TPPI, Echo-AntiEcho**: Different acquisition modes will affect how you construct the multicomplex array. Consult your spectrometer's documentation.

3. **Linear Prediction**: Can be applied before FFT to extend the FID:
   ```julia
   # Apply LP separately to each complex component
   # This would require extracting components, applying LP, and reconstructing
   ```

4. **Baseline Correction**: Usually applied to the magnitude spectrum after all processing.

5. **Integration with NMRTools.jl**: For production NMR processing, see [NMRTools.jl](https://github.com/waudbygroup/NMRTools.jl), which builds on this package with additional NMR-specific functionality.

---

## Performance Considerations

1. **In-place Operations**: `fft!` modifies the array in-place. Copy if you need the original:
   ```julia
   fid_copy = copy(fid)
   fft!(fid_copy, 1, dims=1)
   ```

2. **Memory Layout**: Multicomplex arrays are cache-friendly - all components are stored contiguously.

3. **Type Stability**: Keep types consistent (don't mix Float32/Float64) for best performance.

4. **Preallocate**: Reuse arrays when processing multiple spectra:
   ```julia
   workspace = similar(fid)
   for spectrum in spectra
       workspace .= spectrum
       fft!(workspace, 1, dims=1)
       fft!(workspace, 2, dims=2)
       # ... process workspace ...
   end
   ```

---

## See Also

- **[FFT Operations](@ref fft-operations)**: Detailed FFT documentation
- **[Accessing Components](@ref components)**: Extract parts and components
- **[Examples](@ref Examples)**: Quick code snippets
- **[API Reference](@ref)**: Complete function documentation
- **[NMRTools.jl](https://github.com/waudbygroup/NMRTools.jl)**: Production NMR processing toolkit
