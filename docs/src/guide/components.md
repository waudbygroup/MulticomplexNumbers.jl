# [Accessing Components](@id components)

```@meta
CurrentModule = MulticomplexNumbers
```

This page explains how to access and extract components from multicomplex numbers.


## Real and Imaginary Parts

The `real` and `imag` functions return multicomplex numbers of order N-1:

```@repl components
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

real(z)         # The "real" part with respect to im2
imag(z)         # The "imaginary" part with respect to im2

real(real(z))   # Gets the truly real component
realest(z)      # Get the first (most real) component
```


## Direct Component Access

```@repl components
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

flat(z)          # Get all components as a vector

component(z, 1)  # Get a specific component (1-indexed)
component(z, 2)  # Coefficient of im1
component(z, 3)  # Coefficient of im2
component(z, 4)  # Coefficient of im1*im2
```


## Component Ordering

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


## Order and Properties

```@repl components
using MulticomplexNumbers

order(1.0 + im1)        # Get the order
order(1.0 + im1 + im2)

isreal(5.0 + 0.0*im1)   # Check properties
isreal(1.0 + 2.0*im1)
```


## Matrix Representations

Every multicomplex number can be represented as a real matrix:

```@repl components
using MulticomplexNumbers

z = 1.0 + 2.0*im1
M = matrep(z)

size(M)     # Matrix properties

M[:, 1]     # The first column contains the components
flat(z)

w = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2
size(matrep(w))     # Bicomplex gives 4×4
```

The matrix representation preserves algebra:

```@repl components
w = 3.0 + 4.0*im1

flat(z * w)                     # Direct multiplication
(matrep(z) * matrep(w))[:, 1]   # Multiplication via matrices
```


## Complex Array Views

The `ascomplex` function views multicomplex arrays as complex arrays. This is particularly useful for extracting specific complex sub-components in NMR applications.

```@repl components
using MulticomplexNumbers

data = [Multicomplex(1.0, 2.0, 3.0, 4.0) for i in 1:3]  # Create a bicomplex array

c1 = ascomplex(data, 1)  # View as complex along im1 (extracts real + im1*i components)
size(c1)  # (2, 3) - the 2 is from 2^(N-1)

c2 = ascomplex(data, 2)  # View as complex along im2 (extracts real + im2*i components)
size(c2)

c_default = ascomplex(data)  # Default is highest order, same as c2
```