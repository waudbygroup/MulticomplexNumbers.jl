# Getting Started

This guide will help you get up and running with MulticomplexNumbers.jl quickly.

## Installation

MulticomplexNumbers.jl is a registered Julia package. Install it using Julia's package manager:

```julia
using Pkg
Pkg.add("MulticomplexNumbers")
```

Or in the package REPL (press `]`):
```
pkg> add MulticomplexNumbers
```

### Optional: FFT Support

For Fourier transform capabilities, also install FFTW:
```julia
Pkg.add("FFTW")
```

## Basic Usage

```@repl getting-started
using MulticomplexNumbers
```

### Creating Multicomplex Numbers

The package provides imaginary units `im1` through `im6` for constructing multicomplex numbers:

```@repl getting-started
# Order 1: Complex numbers (2 components)
z1 = 1.0 + 2.0*im1

# Order 2: Bicomplex numbers (4 components)
z2 = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# Order 3: Tricomplex numbers (8 components)
z3 = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im3
```

You can also use direct constructors:

```@repl getting-started
# From real components
Multicomplex(1.0, 2.0)           # = 1 + 2*im1
Multicomplex(1.0, 2.0, 3.0, 4.0) # = 1 + 2*im1 + 3*im2 + 4*im1*im2

# From Complex numbers
Multicomplex(1.0 + 2.0im)        # = 1 + 2*im1
```

### Key Properties

Each imaginary unit squares to -1:
```@repl getting-started
im1 * im1
im2 * im2
```

Unlike quaternions, multicomplex units **commute**:
```@repl getting-started
im1 * im2 == im2 * im1
```

However, the product of different units is hyperbolic (squares to +1):
```@repl getting-started
(im1 * im2)^2
```

### Arithmetic Operations

All standard operations are supported:

```@repl getting-started
a = 1.0 + 2.0*im1 + 3.0*im2
b = 2.0 + 1.0*im1 - 1.0*im2

# Addition and subtraction
a + b
a - b

# Multiplication and division
a * b
a / b

# Powers and transcendental functions
a^2
exp(a)
log(a)
sqrt(a)
```

### Extracting Components

```@repl getting-started
z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# Real and imaginary parts (returns multicomplex of order N-1)
real(z)
imag(z)

# Get all components as a vector
flat(z)

# Get the "most real" component (first component)
realest(z)

# Get a specific component by index
component(z, 2)  # The im1 coefficient
```

## Quick Example: Numerical Differentiation

The most powerful application of multicomplex numbers is computing derivatives with machine precision:

```@repl getting-started
# Compute the derivative of f(x) = sin(x) at x = π/4
f(x) = sin(x)

# Use an extremely small step size (impossible with finite differences!)
h = 1e-100
x = π/4 + h*im1

# Evaluate the function
result = f(x)

# Extract the derivative from the imaginary part
derivative = imag(result) / h

# Compare to exact value: cos(π/4) = √2/2 ≈ 0.7071...
cos(π/4)
```

The result matches the exact derivative to machine precision!

## Next Steps

- **[Background](@ref)**: Learn the mathematics behind multicomplex numbers
- **User Guide**: Detailed exploration of all features:
  - [Creating Numbers](@ref creating)
  - [Arithmetic Operations](@ref arithmetic)
  - [Accessing Components](@ref components)
  - [FFT Operations](@ref fft-operations)
- **Applications**: Practical examples:
  - [NMR Spectroscopy](@ref nmr-applications)
  - [Numerical Differentiation](@ref differentiation-applications)
- **[API Reference](@ref)**: Complete function documentation
