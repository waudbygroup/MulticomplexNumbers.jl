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
z1 = 1.0 + 2.0*im1  # Order 1: Complex numbers (2 components)

z2 = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2  # Order 2: Bicomplex numbers (4 components)

z3 = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im3  # Order 3: Tricomplex numbers (8 components)
```

You can also use direct constructors:

```@repl getting-started
Multicomplex(1.0, 2.0)           # From real components
Multicomplex(1.0, 2.0, 3.0, 4.0)

Multicomplex(1.0 + 2.0im)        # From complex numbers
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

However, the product of different units is **hyperbolic** (squares to +1):
```@repl getting-started
(im1 * im2)^2
```

### Arithmetic Operations

All standard operations are supported:

```@repl getting-started
a = 1.0 + 2.0*im1 + 3.0*im2
b = 2.0 + 1.0*im1 - 1.0*im2

a + b     # Addition
a - b     # Subtraction
a * b     # Multiplication
a / b     # Division
a^2       # Powers
sqrt(a)   # Square root
exp(a)    # Exponentials
log(a)    # Logarithms
```

### Extracting Components

```@repl getting-started
z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

real(z)   # Real / imaginary parts (returns multicomplex of order N-1)
imag(z)

flat(z)   # Get all components as a vector

realest(z)  # Get the "most real" component (first component)

component(z, 2)  # Get a specific component by index
```

## Quick Example: Numerical Differentiation

One application of multicomplex numbers is computing derivatives with machine precision:

```@repl getting-started
f(x) = sin(x)   # Compute the derivative of f(x) = sin(x) at x = π/4

h = 1e-100      # Use an extremely small step size
x = π/4 + h*im1

result = f(x)   # Evaluate the function

derivative = imag(result) / h   # Extract the derivative from the imaginary part

cos(π/4)    # Compare to exact value: cos(π/4) = √2/2 ≈ 0.7071...
```

The result matches the exact derivative to machine precision!