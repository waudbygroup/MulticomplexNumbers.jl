# API Reference

```@meta
CurrentModule = MulticomplexNumbers
```

This page provides complete documentation for all exported types and functions.


## Types

```@docs
Multicomplex
```

## Constants

### Imaginary Units

The package exports imaginary units `im1` through `im6` for conveniently constructing multicomplex numbers.

```@docs
im1
```

Additional units `im2`, `im3`, `im4`, `im5`, and `im6` follow the same pattern for higher orders.

### Programmatic Construction

```@docs
imN
```



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


## Representations
### Matrix Representation

```@docs
matrep
```

### Complex View

```@docs
ascomplex
```


## Standard Functions


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

### Trigonometric Functions

- `sin(m)`: Sine
- `cos(m)`: Cosine
- `tan(m)`: Tangent
- `cot(m)`: Cotangent
- `sec(m)`: Secant
- `csc(m)`: Cosecant

### Hyperbolic Functions

- `sinh(m)`: Hyperbolic sine
- `cosh(m)`: Hyperbolic cosine
- `tanh(m)`: Hyperbolic tangent
- `coth(m)`: Hyperbolic cotangent
- `sech(m)`: Hyperbolic secant
- `csch(m)`: Hyperbolic cosecant

### Inverse Trigonometric Functions

- `asin(m)`: Arcsine
- `acos(m)`: Arccosine
- `atan(m)`: Arctangent
- `acot(m)`: Arccotangent
- `asec(m)`: Arcsecant
- `acsc(m)`: Arccosecant

### Inverse Hyperbolic Functions

- `asinh(m)`: Inverse hyperbolic sine
- `acosh(m)`: Inverse hyperbolic cosine
- `atanh(m)`: Inverse hyperbolic tangent
- `acoth(m)`: Inverse hyperbolic cotangent
- `asech(m)`: Inverse hyperbolic secant
- `acsch(m)`: Inverse hyperbolic cosecant

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


## Multicomplex-Specific Operations

### Fold Operator
```@docs
fold
```

The fold operator multiplies a multicomplex number by its conjugate, reducing the order by one.

### Abient Numbers

```@docs
isabient
```

A multicomplex number is "abient" if it becomes zero after sufficiently many foldings. This property arises because multicomplex numbers are not a composition algebra.

## Random Number Generation

The package supports random number generation for multicomplex numbers:


```julia
rand([rng], Multicomplex{T,N,C})            # single random number
rand([rng], Multicomplex{T,N,C}, dims...)    # array of random numbers
```

Generates random multicomplex numbers with components drawn from the default distribution for type `T`.

**Example:**
```julia
julia> rand(Multicomplex{Float64,1,2})
0.234 + 0.891*im1

julia> rand(Multicomplex{Float64,2,4}, 3)    # vector of 3 random bicomplex numbers
3-element Vector{Multicomplex{Float64, 2, 4}}: ...

julia> rand(Multicomplex{Float64,1,2}, 2, 3) # 2×3 matrix
2×3 Matrix{Multicomplex{Float64, 1, 2}}: ...
```

### Normal Distribution

```julia
randn([rng], Multicomplex{T,N,C})            # single random number
randn([rng], Multicomplex{T,N,C}, dims...)   # array of random numbers
```

Generates random multicomplex numbers with components drawn from a standard normal distribution. Only available for floating-point types.

**Example:**
```julia
julia> randn(Multicomplex{Float64,1,2})
-0.543 + 1.234*im1

julia> randn(Multicomplex{Float64,2,4}, 5)   # vector of 5 random bicomplex numbers
5-element Vector{Multicomplex{Float64, 2, 4}}: ...
```



When FFTW.jl is loaded, the following functions become available. All require an explicit `unit` argument — an integer or a multicomplex imaginary constant (e.g. `im2`).

### In-place Transforms
```julia
fft!(A, unit [, dims])     # Forward FFT
ifft!(A, unit [, dims])    # Inverse FFT (normalized)
bfft!(A, unit [, dims])    # Backward FFT (unnormalized inverse)
```

### Allocating Transforms

```julia
fft(A, unit [, dims])      # Forward FFT (returns copy)
ifft(A, unit [, dims])     # Inverse FFT (returns copy)
bfft(A, unit [, dims])     # Backward FFT (returns copy)
```

**Arguments:**
- `A`: Array of multicomplex numbers
- `unit`: Which imaginary unit to use — an integer (1 for `im1`, 2 for `im2`, etc.) or a multicomplex constant (`im1`, `im2`, etc.)
- `dims`: Dimensions along which to compute the FFT (optional, defaults to all)

**Supported orders:** N = 1, 2, 3, 4

**Example:**
```julia
using MulticomplexNumbers
using FFTW

data = rand(Multicomplex{Float64,2,4}, 64)

# In-place forward and inverse
fft!(data, 1)
ifft!(data, 1)

# Using multicomplex constants
fft!(data, im1)

# Allocating (preserves original)
spectrum = fft(data, 1)
```


## Index
```@index
Pages = ["api.md"]
```
