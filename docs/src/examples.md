# [Examples](@id Examples)

```@meta
CurrentModule = MulticomplexNumbers
```

This page provides quick examples demonstrating key features of MulticomplexNumbers.jl.

---

## Basic Multicomplex Algebra

### Imaginary Units

Multicomplex units `im1` through `im6` are defined. Each squares to -1:

```jldoctest examples
julia> using MulticomplexNumbers

julia> im1
0 + 1*im1

julia> im1 * im1
-1 + 0*im1

julia> im2
(0 + 0*im1) + (1 + 0*im1)*im2

julia> im2 * im2
(-1 + 0*im1) + (0 + 0*im1)*im2
```

Products of different units commute and create new basis elements:

```jldoctest examples
julia> im1 * im2
(0 + 0*im1) + (0 + 1*im1)*im2

julia> im1 * im2 == im2 * im1
true
```

Interestingly, products of distinct units square to +1 (hyperbolic units):

```jldoctest examples
julia> (im1 * im2)^2
(1 + 0*im1) + (0 + 0*im1)*im2
```

### Creating Numbers

```jldoctest examples
julia> z = 1.0 + 2.0*im1  # Complex (order 1)
1.0 + 2.0*im1

julia> w = 3.0 + 4.0*im1 + 5.0*im2 + 6.0*im1*im2  # Bicomplex (order 2)
(3.0 + 4.0*im1) + (5.0 + 6.0*im1)*im2

julia> Multicomplex(1.0, 2.0)  # From components
1.0 + 2.0*im1

julia> Multicomplex(1+2im)  # From Complex
1 + 2*im1
```

---

## Arithmetic

### Basic Operations

```jldoctest examples
julia> a = 1.0 + 2.0*im1;

julia> b = 3.0 + 4.0*im1;

julia> a + b
4.0 + 6.0*im1

julia> a * b
-5.0 + 10.0*im1

julia> a / b
0.44 + 0.08*im1
```

### Powers and Transcendentals

```jldoctest examples
julia> z = 1.0 + 1.0*im1;

julia> z^2
0.0 + 2.0*im1

julia> exp(z)
1.4686939399158851 + 2.2873552871788423*im1

julia> log(z)
0.34657359027997264 + 0.7853981633974483*im1

julia> sqrt(z)
1.0986841134678098 + 0.45508986056222733*im1
```

---

## Component Access

### Extracting Parts

```jldoctest examples
julia> z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2;

julia> real(z)  # Order N-1 real part
1.0 + 2.0*im1

julia> imag(z)  # Order N-1 imaginary part
3.0 + 4.0*im1

julia> flat(z)  # All components
4-element StaticArraysCore.SVector{4, Float64} with indices SOneTo(4):
 1.0
 2.0
 3.0
 4.0

julia> realest(z)  # First component
1.0

julia> component(z, 4)  # Specific component (im1*im2 coefficient)
4.0
```

### Order and Properties

```jldoctest examples
julia> order(1.0 + im1)
1

julia> order(1.0 + im1 + im2)
2

julia> isreal(5.0 + 0.0*im1)
true

julia> isreal(1.0 + 2.0*im1)
false
```

---

## Matrix Representation

Every multicomplex number has an equivalent matrix:

```jldoctest examples
julia> z = 1.0 + 2.0*im1;

julia> M = matrep(z)
2×2 StaticArraysCore.SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 1.0  -2.0
 2.0   1.0
```

The matrix representation preserves algebra:

```jldoctest examples
julia> w = 3.0 + 4.0*im1;

julia> flat(z * w)
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 -5.0
 10.0

julia> (matrep(z) * matrep(w))[:, 1]
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 -5.0
 10.0
```

---

## Numerical Differentiation

### First Derivative

```jldoctest examples
julia> f(x) = x^3;

julia> x = 2.0;

julia> h = 1e-100;

julia> derivative = imag(f(x + h*im1)) / h
12.0

julia> # Exact: f'(x) = 3x² = 12
```

### Second Derivative

```jldoctest examples
julia> f(x) = x^4;

julia> x = 2.0;

julia> h = 1e-50;

julia> result = f(x + h*im1 + h*im2);

julia> second_derivative = component(result, 4) / h^2
48.0

julia> # Exact: f''(x) = 12x² = 48
```

---

## Complex Views

Convert multicomplex arrays to complex arrays:

```jldoctest examples
julia> data = [Multicomplex(1.0, 2.0, 3.0, 4.0)];

julia> c = ascomplex(data, 1);

julia> size(c)
(2, 1)

julia> c[1]
1.0 + 2.0im

julia> c[2]
3.0 + 4.0im
```

---

## Type Promotions

Different orders automatically promote:

```jldoctest examples
julia> a = 1.0 + im1;  # Order 1

julia> b = 2.0 + im2;  # Order 2

julia> c = a + b;  # Promoted to order 2

julia> order(c)
2
```

Scalars and Complex promote to multicomplex:

```jldoctest examples
julia> z = 1.0 + im1;

julia> z + 5.0
6.0 + 1.0*im1

julia> z + (2.0 + 3.0im)
3.0 + 4.0*im1
```

---

## FFT (requires FFTW)

```julia
using MulticomplexNumbers
using FFTW

# Create bicomplex signal
signal = [Multicomplex(cos(2π*k/64), sin(2π*k/64), 0.0, 0.0) for k in 0:63]

# FFT along im1
fft!(signal, 1)

# Peak at frequency 1
abs(realest(signal[2]))  # Should be large
```

---

## High Precision

```julia
julia> z = big"1.0" + big"2.0" * im1;

julia> typeof(z)
Multicomplex{BigFloat, 1, 2}

julia> exp(z)
-1.131204383756813... + 2.471726672004818...*im1
```

---

## Quick Reference

| Task | Code |
|------|------|
| Create order 1 | `1.0 + 2.0*im1` |
| Create order 2 | `1.0 + 2.0*im1 + 3.0*im2` |
| Get components | `flat(z)` |
| Get order | `order(z)` |
| 1st derivative | `imag(f(x + h*im1)) / h` |
| 2nd derivative | `component(f(x + h*im1 + h*im2), 4) / h²` |
| Matrix form | `matrep(z)` |
| As complex | `ascomplex(data, unit)` |
