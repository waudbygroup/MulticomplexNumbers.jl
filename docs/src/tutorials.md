# [Tutorials](@id tutorials)

This section provides practical tutorials for common use cases of multicomplex numbers.

```@meta
CurrentModule = MulticomplexNumbers
```

## Tutorial 1: Numerical Differentiation

The multicomplex step method is the primary application of multicomplex numbers. It extends the complex step derivative to compute higher-order derivatives with machine precision, avoiding the subtractive cancellation errors that plague finite difference methods.

### The Problem with Finite Differences

Traditional numerical differentiation uses finite differences:

```math
f'(x) \approx \frac{f(x+h) - f(x)}{h}
```

This suffers from a fundamental trade-off:
- **Large h**: Truncation error dominates (formula is only approximate)
- **Small h**: Roundoff error dominates (subtractive cancellation in `f(x+h) - f(x)`)

The optimal `h` is typically around `√ε ≈ 10⁻⁸` for Float64, limiting accuracy to about 8 digits.

### The Complex Step Method

For an analytic function `f`, we have the Taylor expansion:

```math
f(x + ih) = f(x) + ih \cdot f'(x) - \frac{h^2}{2}f''(x) - \frac{ih^3}{6}f'''(x) + ...
```

Taking the imaginary part:

```math
\text{Im}[f(x + ih)] = h \cdot f'(x) + O(h^3)
```

Thus:

```math
f'(x) \approx \frac{\text{Im}[f(x + ih)]}{h}
```

Since we're not computing a difference, there's **no subtractive cancellation**! We can use arbitrarily small `h`, giving nearly machine-precision derivatives.

### First Derivatives with `im1`

```@repl tutorials
using MulticomplexNumbers

# Define a test function
f(x) = sin(x) * exp(-x^2)

# Compute derivative at x = 0.5
x = 0.5
h = 1e-100  # Extremely small step - impossible with finite differences!

# Evaluate f at x + h*im1
result = f(x + h*im1)

# Extract derivative from imaginary part
f_prime = imag(result) / h

# Compute exact derivative for comparison: f'(x) = cos(x)*exp(-x²) - 2x*sin(x)*exp(-x²)
exact = cos(x)*exp(-x^2) - 2x*sin(x)*exp(-x^2)

# Error
abs(f_prime - exact)
```

The error is at the level of machine epsilon!

### Second Derivatives with `im1` and `im2`

For second derivatives, we use bicomplex numbers (`im1` and `im2`):

```math
f''(x) \approx \frac{\text{component}_{i_1 i_2}[f(x + h \cdot i_1 + h \cdot i_2)]}{h^2}
```

```@repl tutorials
using MulticomplexNumbers

# Compute second derivative of f(x) = x^4 at x = 2
# Exact answer: f''(x) = 12x² = 48

x = 2.0
h = 1e-50

# Evaluate at x + h*im1 + h*im2
result = (x + h*im1 + h*im2)^4

# The i₁i₂ component is at index 4: [1, i₁, i₂, i₁i₂]
second_derivative = component(result, 4) / h^2

# Exact value
exact = 12 * x^2

abs(second_derivative - exact)
```

### Third Derivatives with `im1`, `im2`, and `im3`

```@repl tutorials
using MulticomplexNumbers

# Compute third derivative of f(x) = x^5 at x = 1
# Exact answer: f'''(x) = 60x² = 60

x = 1.0
h = 1e-30

# Evaluate at x + h*im1 + h*im2 + h*im3
result = (x + h*im1 + h*im2 + h*im3)^5

# The i₁i₂i₃ component is at index 8
third_derivative = component(result, 8) / h^3

exact = 60.0
abs(third_derivative - exact)
```

### Practical Tips for Numerical Differentiation

1. **Step size**: Use `h ≈ 10⁻¹⁰⁰/n` for the n-th derivative. Smaller is generally better.

2. **Component indices**: For an n-th derivative, extract component `2^n` from an order-n multicomplex number.

3. **Computational cost**: Storage and computation scale as `2^n` for the n-th derivative. For orders above ~10, this becomes expensive.

4. **Function requirements**: Your function must accept multicomplex arguments and use operations that work on them (`+`, `-`, `*`, `/`, `^`, `exp`, `log`, `sqrt`, etc.).

### Derivative Helper Function

Here's a convenient wrapper function:

```julia
using MulticomplexNumbers

"""
    derivative(f, x, n=1; h=1e-50)

Compute the n-th derivative of f at x using multicomplex step method.
"""
function derivative(f, x::Real, n::Int=1; h::Real=1e-50)
    if n < 1
        return f(x)
    elseif n == 1
        result = f(x + h*im1)
        return imag(result) / h
    elseif n == 2
        result = f(x + h*im1 + h*im2)
        return component(result, 4) / h^2
    elseif n == 3
        result = f(x + h*im1 + h*im2 + h*im3)
        return component(result, 8) / h^3
    elseif n == 4
        result = f(x + h*im1 + h*im2 + h*im3 + h*im4)
        return component(result, 16) / h^4
    else
        error("Derivatives of order > 4 not implemented in this example")
    end
end

# Test it
f(x) = sin(x)
x = 1.0

derivative(f, x, 1)  # f'(1) = cos(1) ≈ 0.5403
derivative(f, x, 2)  # f''(1) = -sin(1) ≈ -0.8415
derivative(f, x, 3)  # f'''(1) = -cos(1) ≈ -0.5403
derivative(f, x, 4)  # f''''(1) = sin(1) ≈ 0.8415
```

---

## Tutorial 2: Multi-dimensional NMR Signal Processing

Multicomplex numbers provide a natural representation for multi-dimensional NMR data, where each indirect dimension has its own imaginary component. This tutorial shows how to use multicomplex FFT for processing NMR-like signals.

### Background: NMR and Quadrature Detection

In NMR spectroscopy:
- Each time-domain dimension uses **quadrature detection**, recording both real and imaginary components
- A 2D NMR experiment produces data that is complex in both dimensions
- This is naturally represented as bicomplex numbers: `z = a + b*i₁ + c*i₂ + d*i₁i₂`

Where:
- `i₁` represents the imaginary component of the direct (acquisition) dimension
- `i₂` represents the imaginary component of the indirect dimension

### Setting Up FFT Support

```julia
using MulticomplexNumbers
using FFTW  # This loads the FFT extension
```

### Creating Test Signals

```julia
using MulticomplexNumbers
using FFTW

# Create a 2D time-domain signal (bicomplex)
# This simulates a signal with frequencies in both dimensions

function create_2d_signal(n1, n2, freq1, freq2, decay1, decay2)
    t1 = range(0, 1, length=n1)  # Direct dimension time
    t2 = range(0, 1, length=n2)  # Indirect dimension time

    signal = Array{Multicomplex{Float64, 2, 4}}(undef, n1, n2)

    for (j, τ2) in enumerate(t2), (i, τ1) in enumerate(t1)
        # Signal in direct dimension: exp(i*ω₁*t₁) * exp(-t₁/T₂)
        # Signal in indirect dimension: exp(i*ω₂*t₂) * exp(-t₂/T₂)

        # Real/imaginary in dim 1 (im1)
        s1_r = cos(2π * freq1 * τ1) * exp(-τ1 * decay1)
        s1_i = sin(2π * freq1 * τ1) * exp(-τ1 * decay1)

        # Real/imaginary in dim 2 (im2)
        s2_r = cos(2π * freq2 * τ2) * exp(-τ2 * decay2)
        s2_i = sin(2π * freq2 * τ2) * exp(-τ2 * decay2)

        # Combine as bicomplex: (s1_r + s1_i*im1) * (s2_r + s2_i*im2)
        # = s1_r*s2_r + s1_i*s2_r*im1 + s1_r*s2_i*im2 + s1_i*s2_i*im1*im2
        signal[i, j] = Multicomplex(
            s1_r * s2_r,     # real part
            s1_i * s2_r,     # im1 coefficient
            s1_r * s2_i,     # im2 coefficient
            s1_i * s2_i      # im1*im2 coefficient
        )
    end

    return signal
end

# Create a test signal
n = 64
signal = create_2d_signal(n, n, 5.0, 10.0, 2.0, 3.0)
# size(signal) => (64, 64)
# eltype(signal) => Multicomplex{Float64, 2, 4}
```

### FFT Along Specific Dimensions

The `fft!` function performs in-place FFT along a specific imaginary unit:

```julia
using MulticomplexNumbers
using FFTW

# Create a simple bicomplex array for demonstration
data = [Multicomplex(1.0, 0.0, 0.0, 0.0) for i in 1:8]

# FFT along im1 (direct dimension)
fft!(data, 1)
# After FFT along im1: [8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Note: The original data is modified in-place
```

### Complete 2D Transform Example

For a 2D NMR spectrum, you typically:
1. FFT along the direct dimension (im1)
2. FFT along the indirect dimension (im2)

```julia
using MulticomplexNumbers
using FFTW

# Create a 2D bicomplex signal
function make_test_fid()
    n = 32
    data = Array{Multicomplex{Float64, 2, 4}}(undef, n, n)

    for j in 1:n, i in 1:n
        # Simple oscillating signal
        t1, t2 = (i-1)/n, (j-1)/n
        ω1, ω2 = 3.0, 7.0  # Frequencies

        r1, i1 = cos(2π*ω1*t1), sin(2π*ω1*t1)
        r2, i2 = cos(2π*ω2*t2), sin(2π*ω2*t2)

        data[i,j] = Multicomplex(r1*r2, i1*r2, r1*i2, i1*i2)
    end
    return data
end

fid = make_test_fid()

# Transform both dimensions
fid_copy = copy(fid)
fft!(fid_copy, 1)  # Direct dimension
fft!(fid_copy, 2)  # Indirect dimension

# Peak at position (4, 8) contains the signal
abs(realest(fid_copy[4, 8]))  # Should be large
```

### Working with Higher Dimensions

For 3D NMR data, use tricomplex numbers with `im1`, `im2`, `im3`:

```julia
using MulticomplexNumbers
using FFTW

# Create 3D data with tricomplex numbers
data_3d = Array{Multicomplex{Float64, 3, 8}}(undef, 16, 16, 16)

# Fill with test signal...
# (Each element has 8 components: 1, im1, im2, im1*im2, im3, im1*im3, im2*im3, im1*im2*im3)

# FFT each dimension
fft!(data_3d, 1)  # Direct dimension (im1)
fft!(data_3d, 2)  # First indirect (im2)
fft!(data_3d, 3)  # Second indirect (im3)
```

### Viewing Multicomplex Data as Complex

The `ascomplex` function provides views of multicomplex arrays as complex arrays:

```@repl tutorials
using MulticomplexNumbers

# Create bicomplex data
data = [Multicomplex(1.0, 2.0, 3.0, 4.0) for _ in 1:4]

# View as complex along im1
view_im1 = ascomplex(data, 1)
size(view_im1)

# View as complex along im2
view_im2 = ascomplex(data, 2)
size(view_im2)
```

### Tips for NMR Applications

1. **Choose the right order**: Use bicomplex (N=2) for 2D, tricomplex (N=3) for 3D, etc.

2. **Component ordering**: Components follow a binary pattern:
   - Order 2: `[1, i₁, i₂, i₁i₂]` (indices 1-4)
   - Order 3: `[1, i₁, i₂, i₁i₂, i₃, i₁i₃, i₂i₃, i₁i₂i₃]` (indices 1-8)

3. **Integration with NMRTools.jl**: See [NMRTools.jl](https://github.com/waudbygroup/NMRTools.jl) for comprehensive NMR data handling that builds on this package.

4. **Memory efficiency**: Multicomplex arrays store components contiguously, making them cache-friendly for numerical operations.

---

## Tutorial 3: Matrix Representations

Every multicomplex number has an equivalent matrix representation that preserves algebraic operations. This is useful for understanding the structure and for implementing functions like `exp`, `log`, and `sqrt`.

### The Matrix Representation

For complex numbers, the 2×2 matrix representation is:

```math
a + bi \mapsto \begin{pmatrix} a & -b \\ b & a \end{pmatrix}
```

For bicomplex numbers, this extends to 4×4 matrices, and so on.

```@repl tutorials
using MulticomplexNumbers

# Complex number matrix representation
z = 3.0 + 4.0*im1
M_z = matrep(z)

# Verify: multiplication via matrices
w = 1.0 + 2.0*im1
M_w = matrep(w)

# Matrix product corresponds to multicomplex product
z * w
Multicomplex{1}(M_z * M_w * [1, 0])  # First column of product
```

### Using Matrix Representations

```@repl tutorials
using MulticomplexNumbers

# Bicomplex number
z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2
M = matrep(z)

# The matrix is 4×4 for bicomplex
size(M)

# First column contains the components
M[:, 1]
flat(z)

# Matrix functions give multicomplex functions
using LinearAlgebra
exp_M = exp(M)
exp_z = exp(z)

# Compare first columns
exp_M[:, 1]
flat(exp_z)
```

This matrix approach is how `exp`, `log`, `sqrt`, and `^` are implemented internally.

---

## Summary

| Application | Order | Units Used | Key Functions |
|-------------|-------|------------|---------------|
| 1st derivative | 1 | `im1` | `imag(f(x + h*im1)) / h` |
| 2nd derivative | 2 | `im1`, `im2` | `component(result, 4) / h²` |
| 3rd derivative | 3 | `im1`, `im2`, `im3` | `component(result, 8) / h³` |
| 2D NMR | 2 | `im1`, `im2` | `fft!(data, 1)`, `fft!(data, 2)` |
| 3D NMR | 3 | `im1`, `im2`, `im3` | `fft!(data, 1)`, `fft!(data, 2)`, `fft!(data, 3)` |

## References

- [Complex Step Approximation (Nick Higham)](https://nhigham.com/2020/10/06/what-is-the-complex-step-approximation/)
- Lantoine, G., Russell, R. P., Dargent, T. "Using Multicomplex Variables for Automatic Computation of High-order Derivatives." ACM TOMS, 2012.
- [NIST Multicomplex Report](https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf)
