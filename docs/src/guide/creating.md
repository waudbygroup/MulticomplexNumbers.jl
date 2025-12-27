# Creating Multicomplex Numbers

```@meta
CurrentModule = MulticomplexNumbers
```

This page explains how to create and work with multicomplex numbers.

---

## The Multicomplex Type

The core type is [`Multicomplex{T,N,C}`](@ref) where:
- `T`: The base scalar type (must be `<:Real`, e.g., `Float64`, `BigFloat`)
- `N`: The order (number of imaginary units)
- `C`: The number of components (always `C = 2^N`)

```@repl creating
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

## Using Imaginary Units

The package exports [`im1`](@ref) through [`im6`](@ref):

```@repl creating
using MulticomplexNumbers

im1  # First imaginary unit
im2  # Second imaginary unit
im3  # Third imaginary unit

# Build multicomplex numbers by combination
z = 3.0 + 4.0*im1                           # Order 1
w = 1.0 + 2.0*im1 + 3.0*im2                 # Order 2 (auto-promoted)
v = 1.0 + im1 + im2 + im3 + im1*im2*im3    # Order 3
```

---

## Direct Constructors

```@repl creating
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

```@repl creating
using MulticomplexNumbers
using StaticArrays

# Construct with explicit type parameters
z = Multicomplex{Float64,2,4}(1.0)  # Bicomplex 1.0

# From StaticArrays SVector
v = SVector(1.0, 2.0, 3.0, 4.0)
Multicomplex{2}(v)
```

---

## Type Conversions and Promotion

### Converting to Multicomplex

```@repl creating
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

```@repl creating
using MulticomplexNumbers

z = Multicomplex(5.0)  # Order 0

# Convert to real if purely real
Float64(z)

# This would error for non-real multicomplex:
# Float64(1.0 + 2.0*im1)  # InexactError
```

---

## Working with Higher Precision

MulticomplexNumbers supports any `Real` scalar type:

```@repl creating
using MulticomplexNumbers

# BigFloat for high precision
z = big"1.0" + big"2.0"*im1
typeof(z)

# Compute with high precision
exp(z)
```

---

## See Also

- **[Arithmetic Operations](@ref arithmetic)**: Mathematical operations on multicomplex numbers
- **[Accessing Components](@ref components)**: Extract parts and components
- **[API Reference](@ref)**: Complete function documentation
