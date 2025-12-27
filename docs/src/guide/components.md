# [Accessing Components](@id components)

```@meta
CurrentModule = MulticomplexNumbers
```

This page explains how to access and extract components from multicomplex numbers.

---

## Real and Imaginary Parts

The `real` and `imag` functions return multicomplex numbers of order N-1:

```@repl components
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# For order 2: real/imag return order 1
real(z)  # The "real" part with respect to im2
imag(z)  # The "imaginary" part with respect to im2

# Apply recursively
real(real(z))  # Gets the truly real component
```

---

## Direct Component Access

```@repl components
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# Get all components as a vector
flat(z)

# Get a specific component (1-indexed)
component(z, 1)  # Coefficient of 1
component(z, 2)  # Coefficient of im1
component(z, 3)  # Coefficient of im2
component(z, 4)  # Coefficient of im1*im2

# Get the first (most real) component
realest(z)
```

---

## Component Ordering

Components follow a binary pattern based on which imaginary units are present:

| Order | Index | Basis Element |
|-------|-------|---------------|
| 1 | 1 | 1 |
| 1 | 2 | i₁ |
| 2 | 1 | 1 |
| 2 | 2 | i₁ |
| 2 | 3 | i₂ |
| 2 | 4 | i₁i₂ |
| 3 | 1 | 1 |
| 3 | 2 | i₁ |
| 3 | 3 | i₂ |
| 3 | 4 | i₁i₂ |
| 3 | 5 | i₃ |
| 3 | 6 | i₁i₃ |
| 3 | 7 | i₂i₃ |
| 3 | 8 | i₁i₂i₃ |

The pattern: index `k` has the units corresponding to bits set in `k-1` (0-indexed).

---

## Order and Properties

```@repl components
using MulticomplexNumbers

# Get the order
order(1.0 + im1)
order(1.0 + im1 + im2)

# Check properties
isreal(5.0 + 0.0*im1)
isreal(1.0 + 2.0*im1)
```

---

## Matrix Representations

Every multicomplex number can be represented as a real matrix:

```@repl components
using MulticomplexNumbers

z = 1.0 + 2.0*im1
M = matrep(z)

# Matrix properties
size(M)

# The first column contains the components
M[:, 1]
flat(z)

# Bicomplex gives 4×4
w = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2
size(matrep(w))
```

The matrix representation preserves algebra:

```@repl components
# Multiplication via matrices
w = 3.0 + 4.0*im1

flat(z * w)
(matrep(z) * matrep(w))[:, 1]
```

---

## Complex Array Views

The `ascomplex` function views multicomplex arrays as complex arrays. This is particularly useful for extracting specific complex sub-components in NMR applications.

```@repl components
using MulticomplexNumbers

# Create a bicomplex array
data = [Multicomplex(1.0, 2.0, 3.0, 4.0) for i in 1:3]

# View as complex along im1 (extracts real + im1*i components)
c1 = ascomplex(data, 1)
size(c1)  # (2, 3) - the 2 is from 2^(N-1)

# View as complex along im2 (extracts real + im2*i components)
c2 = ascomplex(data, 2)
size(c2)

# Default is highest order
c_default = ascomplex(data)  # Same as ascomplex(data, 2)
```

### Extracting Specific Complex Components

For NMR data, you often want to extract specific complex sub-spectra:

```@repl components
using MulticomplexNumbers

# Create bicomplex data
z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# Extract the real/im1 complex component
# This gives you a Complex number from the "real" part (with respect to im2)
real_im1_component = real(z)  # Returns 1.0 + 2.0*im1

# Extract the real/im2 complex component
# Re-arrange to view im2 as the complex unit
# For a single number, you can manually extract:
real_part = realest(z)  # 1.0
im2_part = component(z, 3)  # 3.0
# Complex(real_part, im2_part) would give 1.0 + 3.0im

# For arrays, ascomplex does this automatically
data = [Multicomplex(1.0, 2.0, 3.0, 4.0) for _ in 1:4]
complex_view_im1 = ascomplex(data, 1)  # View as complex along im1
complex_view_im2 = ascomplex(data, 2)  # View as complex along im2
```

---

## Performance Tips

1. **Use `flat()` once**: If you need multiple components, call `flat()` once and index the result

2. **Avoid order changes**: Mixing orders triggers promotion, which allocates

3. **Pre-allocate arrays**: For FFT operations, reuse arrays when possible

4. **Consider StaticArrays**: The internal representation uses `SVector`, which is stack-allocated

5. **Profile high orders**: Order N uses 2^N components. Order 6 has 64 components per number!

---

## See Also

- **[Creating Numbers](@ref creating)**: How to create multicomplex numbers
- **[Arithmetic Operations](@ref arithmetic)**: Mathematical operations
- **[FFT Operations](@ref fft-operations)**: Fourier transforms with NMR applications
- **[API Reference](@ref)**: Complete function documentation
