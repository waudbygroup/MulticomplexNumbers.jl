# API Reference

```@meta
CurrentModule = MulticomplexNumbers
```

This page provides complete documentation for all exported types and functions.

---

## Types

```@docs
Multicomplex
```

---

## Constants

### Imaginary Units

The package exports imaginary units `im1` through `im6` for conveniently constructing multicomplex numbers.

```@docs
im1
```

Additional units `im2`, `im3`, `im4`, `im5`, and `im6` follow the same pattern for higher orders.

---

## Component Access

### Order

```@docs
order
```

### Flat Components

```@docs
flat
```

### Individual Components

```@docs
component
```

### Real Component

```@docs
realest
```

---

## Representations

### Matrix Representation

```@docs
matrep
```

### Complex View

```@docs
ascomplex
```

---

## Standard Functions

The following `Base` functions are extended for multicomplex numbers:

### Real and Imaginary Parts

- `real(m::Multicomplex)`: Returns the real part with respect to the highest imaginary unit (order N-1)
- `imag(m::Multicomplex)`: Returns the imaginary part with respect to the highest imaginary unit (order N-1)

### Arithmetic

- `+`, `-`: Addition and subtraction (element-wise on components)
- `*`: Multiplication (via matrix representation)
- `/`: Division (via matrix representation)
- `^`: Power (integer and real exponents)
- `inv(m)`: Multiplicative inverse

### Transcendental Functions

- `exp(m)`: Exponential
- `log(m)`: Natural logarithm
- `sqrt(m)`: Square root

### Conjugation and Norms

- `conj(m)`: Multicomplex conjugation (negates highest imaginary part)
- `abs(m)`: Absolute value (Euclidean norm)
- `abs2(m)`: Squared absolute value

### Comparisons and Properties

- `==`, `isequal`: Equality comparison
- `isreal(m)`: Check if all imaginary components are zero
- `isinteger(m)`: Check if real and integer-valued
- `isfinite(m)`: Check if all components are finite
- `isnan(m)`: Check if any component is NaN
- `isinf(m)`: Check if any component is infinite
- `iszero(m)`: Check if all components are zero
- `isone(m)`: Check if equal to multiplicative identity

### Type Utilities

- `zero(Multicomplex{T,N,C})`: Additive identity
- `one(Multicomplex{T,N,C})`: Multiplicative identity
- `float(Multicomplex{T,N,C})`: Convert to floating-point base type

---

## FFT Support (FFTW Extension)

When FFTW.jl is loaded, the following function becomes available:

### In-place FFT

```julia
fft!(A::AbstractArray{<:Multicomplex}, unit::Integer)
fft!(A::AbstractArray{<:Multicomplex}, unit::Integer, dims)
```

Performs an in-place FFT on a multicomplex array, treating `unit` as the imaginary unit for the complex FFT.

**Arguments:**
- `A`: Array of multicomplex numbers
- `unit`: Which imaginary unit to use (1 for `im1`, 2 for `im2`, etc.)
- `dims`: Dimensions along which to compute the FFT (optional, defaults to all)

**Supported orders:** N = 1, 2, 3, 4

**Example:**
```julia
using MulticomplexNumbers
using FFTW

data = [Multicomplex(rand(4)...) for _ in 1:64]
fft!(data, 1)  # FFT along im1
fft!(data, 2)  # FFT along im2
```

---

## Index

```@index
Pages = ["api.md"]
```
