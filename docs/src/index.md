# MulticomplexNumbers.jl

A [Julia](http://julialang.org) package for representing [multicomplex numbers](https://en.wikipedia.org/wiki/Multicomplex_number) and performing multicomplex algebra.

## What are Multicomplex Numbers?

Multicomplex numbers are a generalization of complex numbers, recursively defined to contain multiple imaginary units (``i_1``, ``i_2``, ``i_3``, ...). Unlike quaternions and Clifford algebras, these imaginary units **commute**: ``i_1 i_2 = i_2 i_1``.

```math
\mathbb{C}_0 = \mathbb{R}, \quad \mathbb{C}_n = \mathbb{C}_{n-1} + i_n \mathbb{C}_{n-1}
```

| Order | Name | Components | Example |
|-------|------|------------|---------|
| 0 | Real | 1 | `5.0` |
| 1 | Complex | 2 | `1 + 2i₁` |
| 2 | Bicomplex | 4 | `1 + 2i₁ + 3i₂ + 4i₁i₂` |
| 3 | Tricomplex | 8 | `1 + 2i₁ + 3i₂ + 4i₃ + ...` |

## Key Applications

### Multi-dimensional NMR Spectroscopy

Bicomplex and higher-order multicomplex numbers naturally represent multi-dimensional NMR data where each dimension has its own imaginary component:

```julia
using MulticomplexNumbers
using FFTW

# 2D NMR: bicomplex numbers (im1 = direct, im2 = indirect)
# Create bicomplex FID where each dimension has real + imaginary components
fid = Array{Multicomplex{Float64, 2, 4}}(undef, 128, 64)
# ... load your NMR data ...

# Transform each dimension with its associated imaginary unit
fft!(fid, 1, dims=1)  # Direct dimension (im1)
fft!(fid, 2, dims=2)  # Indirect dimension (im2)

# Phase correction using appropriate imaginary units
fid .*= exp(im1 * 0.1)  # Correct direct dimension
fid .*= exp(im2 * 0.2)  # Correct indirect dimension
```

### High-Precision Numerical Differentiation

The multicomplex step method computes derivatives with **machine precision**, avoiding the subtractive cancellation errors that limit finite difference methods:

```julia
using MulticomplexNumbers

# First derivative with machine precision
f(x) = sin(x) * exp(-x^2)
h = 1e-100  # Impossibly small for finite differences!
x = 0.5 + h*im1
derivative = imag(f(x)) / h  # Exact to ~15 digits
```

## Quick Start

```julia
using Pkg
Pkg.add("MulticomplexNumbers")

using MulticomplexNumbers

# Create multicomplex numbers
z = 1.0 + 2.0*im1                    # Complex (order 1)
w = 1.0 + 2.0*im1 + 3.0*im2          # Bicomplex (order 2)

# Arithmetic
z * w
exp(w)
sqrt(w)

# Numerical differentiation
h = 1e-100
x = 2.0 + h*im1
f_prime = imag(x^3) / h  # = 12.0 (exact!)
```

## Features

- **Full arithmetic**: `+`, `-`, `*`, `/`, `^`, `exp`, `log`, `sqrt`
- **Imaginary units**: `im1` through `im6` for orders 1-6
- **Type-stable**: Uses StaticArrays for stack-allocated, cache-friendly storage
- **FFT support**: In-place multicomplex FFT via FFTW extension
- **High precision**: Works with `BigFloat` and other `Real` types
- **Julia integration**: Proper `Number` subtype with promotion rules

## Documentation

```@contents
Pages = [
    "getting-started.md",
    "background.md",
    "guide/creating.md",
    "guide/arithmetic.md",
    "guide/components.md",
    "guide/fft.md",
    "applications/nmr.md",
    "applications/differentiation.md",
    "examples.md",
    "api.md"
]
Depth = 2
```

---

## Installation

MulticomplexNumbers.jl is a [registered package](http://pkg.julialang.org):

```julia
using Pkg
Pkg.add("MulticomplexNumbers")
```

For FFT support:
```julia
Pkg.add("FFTW")
```

---

## Related Packages

- [NMRTools.jl](https://github.com/waudbygroup/NMRTools.jl): A library for handling NMR data in Julia (uses this package)
- [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl): Dual numbers for first-order automatic differentiation
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl): Forward mode automatic differentiation
- [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl): Hyper-dual numbers for second derivatives
- [CliffordAlgebras.jl](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl): Clifford and geometric algebras

---

## References

1. **NIST Report**: Bell, I. H., Deiters, U. K. (2021). "Precise Numerical Differentiation of Thermodynamic Functions with Multicomplex Variables." [NIST J. Res. 126:033](https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf)

2. **ACM Algorithm**: Casado, J. M. V., Hewson, R. (2020). "Algorithm 1008: Multicomplex Number Class for Matlab." ACM Trans. Math. Softw. 46:1-26. [doi:10.1145/3378542](http://dx.doi.org/10.1145/3378542)

3. **NIST Implementation**: [github.com/usnistgov/multicomplex](https://github.com/usnistgov/multicomplex) - C++ multicomplex library

4. **Original Method**: Lantoine, G., Russell, R. P., Dargent, T. (2012). "Using Multicomplex Variables for Automatic Computation of High-order Derivatives." ACM Trans. Math. Softw. 38(3):16. [doi:10.1145/2168773.2168774](https://doi.org/10.1145/2168773.2168774)

---

## Authors

- [Chris Waudby](https://waudbylab.org), UCL School of Pharmacy, London (UK)

## Citing

If you find this package useful, please cite:

> Waudby, C. A. (2022). MulticomplexNumbers.jl. [github.com/waudbygroup/MulticomplexNumbers.jl](https://github.com/waudbygroup/MulticomplexNumbers.jl)

## License

MulticomplexNumbers.jl is licensed under the MIT license; see [LICENSE](https://github.com/waudbygroup/MulticomplexNumbers.jl/blob/main/LICENSE) for the full license text.
