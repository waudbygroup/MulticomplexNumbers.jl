# [Creating Multicomplex Numbers](@id creating)

```@meta
CurrentModule = MulticomplexNumbers
```

This page explains how to create and work with multicomplex numbers.


## The Multicomplex Type

The core type is [`Multicomplex{T,N,C}`](@ref) where:
- `T`: The base scalar type (must be `<:Real`, e.g., `Float64`, `BigFloat`)
- `N`: The order (number of imaginary units)
- `C`: The number of components (always `C = 2^N`)

```@repl creating
using MulticomplexNumbers

z = 1.0 + 2.0*im1
typeof(z)       # Examine type

w = 1.0 + 2.0*im1 + 3.0*im2
typeof(w)       # Examine type
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


## Using Imaginary Units

The package exports `im1` through `im6` as predefined constants, and [`imN(n)`](@ref) for programmatic construction of any imaginary unit:

```@repl creating
using MulticomplexNumbers

im1  # First imaginary unit
im2  # Second imaginary unit
im3  # Third imaginary unit

z = 3.0 + 4.0*im1                          # Order 1
w = 1.0 + 2.0*im1 + 3.0*im2                # Order 2 (auto-promoted)
v = 1.0 + im1 + im2 + im3 + im1*im2*im3    # Order 3
```

For orders beyond 6, or when constructing units programmatically (e.g. in a loop), use [`imN`](@ref):

```@repl creating
using MulticomplexNumbers

imN(2) == im2  # Same as the predefined constant

for n in 1:3
    # Useful in loops over dimensions
    println("im$n squares to ", imN(n)^2 |> realest)
end
```


## Direct Constructors

```@repl creating
using MulticomplexNumbers

Multicomplex(5.0)                       # Order 0 (scalar)
Multicomplex(1.0, 2.0)                  # Order 1 (from two reals)
Multicomplex(3.0 + 4.0im)               # Order 1 (from Complex)
Multicomplex(1.0, 2.0, 3.0, 4.0)        # Order 2 (from four reals)
Multicomplex(1.0 + 2.0im, 3.0 + 4.0im)  # Order 2 (from two Complex numbers)
```

### Type-Specific Constructors

```@repl creating
using MulticomplexNumbers, StaticArrays

z = Multicomplex{Float64,2,4}(1.0)  # Construct with explicit type parameters

v = SVector(1.0, 2.0, 3.0, 4.0)
Multicomplex{2}(v)   # From StaticArrays SVector
```


## Type Conversions and Promotion

### Converting to Multicomplex

```@repl creating
using MulticomplexNumbers

z = 1.0 + 2.0*im1
z + 3.0              # Reals promote to multicomplex

z + (3.0 + 4.0im)    # Complex promotes too

a = 1.0 + im1        
b = 2.0 + im2        
a + b                # Different orders promote to the higher order
```

### Converting from Multicomplex

```@repl creating
using MulticomplexNumbers

z = Multicomplex(5.0)  # Order 0

Float64(z)  # Convert to real if purely real
```


## Working with Higher Precision

MulticomplexNumbers supports any `Real` scalar type:

```@repl creating
using MulticomplexNumbers

z = big"1.0" + big"2.0"*im1     # BigFloat for high precision
typeof(z)
exp(z)      # Compute with high precision
```


## Random Multicomplex Numbers

Generate random multicomplex numbers using `rand` and `randn`:

```@repl creating
using MulticomplexNumbers

rand(Multicomplex{Float64,1,2})       # Uniform components in [0, 1)
randn(Multicomplex{Float64,2,4})      # Normal components
rand(Multicomplex{Float64,1,2}, 3)    # Array of 3 random order-1 numbers
rand(Multicomplex{Float64,2,4}, 2, 3) # 2×3 matrix of random order-2 numbers
```
