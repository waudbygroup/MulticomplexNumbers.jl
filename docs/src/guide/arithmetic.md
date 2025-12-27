# [Arithmetic Operations](@id arithmetic)

```@meta
CurrentModule = MulticomplexNumbers
```

This page covers all arithmetic operations supported by multicomplex numbers.

---

## Basic Arithmetic

```@repl arithmetic
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

---

## Multiplication

Multicomplex multiplication uses the matrix representation internally:

```@repl arithmetic
using MulticomplexNumbers

a = 1.0 + 2.0*im1
b = 3.0 + 4.0*im1

a * b  # = (1*3 - 2*4) + (1*4 + 2*3)*im1 = -5 + 10*im1

# Verify the algebra
im1 * im1  # = -1
im1 * im2  # = im1*im2 (distinct product)
im2 * im1  # Commutative!
```

---

## Division

Division works except when dividing by zero divisors:

```@repl arithmetic
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

---

## Powers

```@repl arithmetic
using MulticomplexNumbers

z = 1.0 + 0.5*im1

z^2
z^3
z^0.5  # Non-integer powers work too
z^(-1)  # Same as inv(z)
```

---

## Transcendental Functions

```@repl arithmetic
using MulticomplexNumbers

z = 0.5 + 0.25*im1

exp(z)
log(z)
sqrt(z)

# Verify: exp(log(z)) ≈ z
exp(log(z))
```

---

## Conjugation

Conjugation inverts the sign of the highest-order imaginary component:

```@repl arithmetic
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

conj(z)  # Inverts im2 component sign

# For order 1, this matches complex conjugation
w = 1.0 + 2.0*im1
conj(w)
```

---

## Absolute Values and Norms

```@repl arithmetic
using MulticomplexNumbers
using LinearAlgebra

z = 3.0 + 4.0*im1

abs2(z)  # Sum of squares of components
abs(z)   # sqrt(abs2(z))
norm(z)  # Euclidean norm (same as abs for order 1)
```

---

## Testing Properties

```@repl arithmetic
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

## See Also

- **[Creating Numbers](@ref creating)**: How to create multicomplex numbers
- **[Accessing Components](@ref components)**: Extract parts and components
- **[API Reference](@ref)**: Complete function documentation
