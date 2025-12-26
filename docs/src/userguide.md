# [User Guide](@id userguide)

```@meta
CurrentModule = MulticomplexNumbers
```

This guide provides a complete tour of MulticomplexNumbers.jl functionality.

---

## The Multicomplex Type

The core type is `Multicomplex{T,N,C}` where:
- `T`: The base scalar type (must be `<:Real`, e.g., `Float64`, `BigFloat`)
- `N`: The order (number of imaginary units)
- `C`: The number of components (always `C = 2^N`)

```@repl userguide
using MulticomplexNumbers

# Examine types
z = 1.0 + 2.0*im1
typeof(z)

w = 1.0 + 2.0*im1 + 3.0*im2
typeof(w)
```

### Type Hierarchy

```
Number
└── Multicomplex{T,N,C}
    ├── Multicomplex{Float64,0,1}  # Just a real number
    ├── Multicomplex{Float64,1,2}  # Complex-like
    ├── Multicomplex{Float64,2,4}  # Bicomplex
    ├── Multicomplex{Float64,3,8}  # Tricomplex
    └── ...
```

---

## Creating Multicomplex Numbers

### Using Imaginary Units

The package exports `im1` through `im6`:

```@repl userguide
using MulticomplexNumbers

im1  # First imaginary unit
im2  # Second imaginary unit
im3  # Third imaginary unit

# Build multicomplex numbers by combination
z = 3.0 + 4.0*im1                           # Order 1
w = 1.0 + 2.0*im1 + 3.0*im2                 # Order 2 (auto-promoted)
v = 1.0 + im1 + im2 + im3 + im1*im2*im3    # Order 3
```

### Direct Constructors

```@repl userguide
using MulticomplexNumbers

# Order 0 (scalar)
Multicomplex(5.0)

# Order 1 (from two reals)
Multicomplex(1.0, 2.0)  # = 1 + 2*im1

# Order 1 (from Complex)
Multicomplex(3.0 + 4.0im)  # = 3 + 4*im1

# Order 2 (from four reals): 1 + 2*im1 + 3*im2 + 4*im1*im2
Multicomplex(1.0, 2.0, 3.0, 4.0)

# Order 2 (from two Complex numbers)
Multicomplex(1.0 + 2.0im, 3.0 + 4.0im)  # = (1+2*im1) + (3+4*im1)*im2
```

### Type-Specific Constructors

```@repl userguide
using MulticomplexNumbers
using StaticArrays

# Construct with explicit type parameters
z = Multicomplex{Float64,2,4}(1.0)  # Bicomplex 1.0

# From StaticArrays SVector
v = SVector(1.0, 2.0, 3.0, 4.0)
Multicomplex{2}(v)
```

---

## Accessing Components

### Real and Imaginary Parts

The `real()` and `imag()` functions return multicomplex numbers of order N-1:

```@repl userguide
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# For order 2: real/imag return order 1
real(z)  # The "real" part with respect to im2
imag(z)  # The "imaginary" part with respect to im2

# Apply recursively
real(real(z))  # Gets the truly real component
```

### Direct Component Access

```@repl userguide
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

### Component Ordering

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

## Arithmetic Operations

### Basic Arithmetic

```@repl userguide
using MulticomplexNumbers

a = 1.0 + 2.0*im1 + 3.0*im2
b = 2.0 - 1.0*im1 + 0.5*im2

# Addition and subtraction
a + b
a - b
-a

# Scalar operations
3.0 * a
a / 2.0
a + 5.0
```

### Multiplication

Multicomplex multiplication uses the matrix representation internally:

```@repl userguide
using MulticomplexNumbers

a = 1.0 + 2.0*im1
b = 3.0 + 4.0*im1

a * b  # = (1*3 - 2*4) + (1*4 + 2*3)*im1 = -5 + 10*im1

# Verify the algebra
im1 * im1  # = -1
im1 * im2  # = im1*im2 (distinct product)
im2 * im1  # Commutative!
```

### Division

Division works except when dividing by zero divisors:

```@repl userguide
using MulticomplexNumbers

a = 1.0 + 2.0*im1
b = 3.0 + 4.0*im1

a / b
inv(b)  # = 1/b
```

!!! warning "Zero Divisors"
    For N ≥ 2, some non-zero multicomplex numbers cannot be inverted:
    ```julia
    z = 1.0 + 1.0*im1*im2  # This has a zero divisor issue
    # (1 + im1*im2)(1 - im1*im2) = 1 - (im1*im2)² = 1 - 1 = 0
    ```

### Powers

```@repl userguide
using MulticomplexNumbers

z = 1.0 + 0.5*im1

z^2
z^3
z^0.5  # Non-integer powers work too
z^(-1)  # Same as inv(z)
```

### Transcendental Functions

```@repl userguide
using MulticomplexNumbers

z = 0.5 + 0.25*im1

exp(z)
log(z)
sqrt(z)

# Verify: exp(log(z)) ≈ z
exp(log(z))
```

---

## Conjugation and Norms

### Multicomplex Conjugation

Conjugation inverts the sign of the highest-order imaginary component:

```@repl userguide
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

conj(z)  # Inverts im2 component sign

# For order 1, this matches complex conjugation
w = 1.0 + 2.0*im1
conj(w)
```

### Absolute Values and Norms

```@repl userguide
using MulticomplexNumbers
using LinearAlgebra

z = 3.0 + 4.0*im1

abs2(z)  # Sum of squares of components
abs(z)   # sqrt(abs2(z))
norm(z)  # Euclidean norm (same as abs for order 1)
```

---

## Type Conversions and Promotion

### Converting to Multicomplex

```@repl userguide
using MulticomplexNumbers

# Reals promote to multicomplex
z = 1.0 + 2.0*im1
z + 3.0  # 3.0 is promoted to Multicomplex

# Complex promotes too
z + (3.0 + 4.0im)

# Different orders promote to the higher order
a = 1.0 + im1        # Order 1
b = 2.0 + im2        # Order 2
a + b                 # Result is order 2
```

### Converting from Multicomplex

```@repl userguide
using MulticomplexNumbers

z = Multicomplex(5.0)  # Order 0

# Convert to real if purely real
Float64(z)

# This would error for non-real multicomplex:
# Float64(1.0 + 2.0*im1)  # InexactError
```

---

## Testing Properties

```@repl userguide
using MulticomplexNumbers

z = 1.0 + 2.0*im1
w = 5.0 + 0.0*im1

isreal(z)   # false
isreal(w)   # true (imaginary part is zero)

isinteger(w)  # true
isinteger(z)  # false

iszero(z)     # false
iszero(0.0 + 0.0*im1)  # true

isone(1.0 + 0.0*im1)   # true

isfinite(z)  # true
isnan(z)     # false
isinf(z)     # false
```

---

## Matrix Representations

Every multicomplex number can be represented as a real matrix:

```@repl userguide
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

---

## Complex Array Views

The `ascomplex` function views multicomplex arrays as complex arrays:

```@repl userguide
using MulticomplexNumbers

# Create a bicomplex array
data = [Multicomplex(1.0, 2.0, 3.0, 4.0) for i in 1:3]

# View as complex along im1
c1 = ascomplex(data, 1)
size(c1)  # (2, 3) - the 2 is from 2^(N-1)

# View as complex along im2
c2 = ascomplex(data, 2)
size(c2)

# Default is highest order
c_default = ascomplex(data)  # Same as ascomplex(data, 2)
```

---

## FFT Support

When FFTW is loaded, in-place FFT is available:

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

See the [Tutorials](@ref tutorials) for complete NMR examples.

---

## Working with Higher Precision

MulticomplexNumbers supports any `Real` scalar type:

```@repl userguide
using MulticomplexNumbers

# BigFloat for high precision
z = big"1.0" + big"2.0"*im1
typeof(z)

# Compute with high precision
exp(z)
```

---

## Performance Tips

1. **Use Float64**: It's the fastest for most applications

2. **Avoid order changes**: Mixing orders triggers promotion, which allocates

3. **Pre-allocate arrays**: For FFT operations, reuse arrays when possible

4. **Consider StaticArrays**: The internal representation uses `SVector`, which is stack-allocated

5. **Profile high orders**: Order N uses 2^N components. Order 6 has 64 components per number!

---

## Common Patterns

### Pattern 1: Numerical Differentiation

```@repl userguide
using MulticomplexNumbers

f(x) = sin(x) * cos(x)
x = 1.0
h = 1e-100

# First derivative
f_prime = imag(f(x + h*im1)) / h

# Exact: f'(x) = cos(2x)
cos(2x)
```

### Pattern 2: Processing 2D Data

```julia
using MulticomplexNumbers
using FFTW

# Create 2D bicomplex data
data = [Multicomplex(rand(4)...) for i in 1:64, j in 1:64]

# FFT both dimensions
fft!(data, 1)  # Direct
fft!(data, 2)  # Indirect
```

### Pattern 3: Working with Components

```@repl userguide
using MulticomplexNumbers

# Extract the im1*im2 coefficient from a result
z = (1.0 + im1 + im2)^2
component(z, 4)  # The im1*im2 term
```

---

## See Also

- **[Background](@ref)**: Mathematical foundations
- **[Tutorials](@ref tutorials)**: Practical applications
- **[API Reference](@ref)**: Complete function documentation
