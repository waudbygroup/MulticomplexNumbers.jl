# Background

---

## Introduction

[MulticomplexNumbers.jl](https://github.com/waudbygroup/MulticomplexNumbers.jl) provides a Julia implementation
of [multicomplex numbers](https://en.wikipedia.org/wiki/Multicomplex_number). These numbers are a generalisation
of complex numbers, recursively defined to contain multiple imaginary numbers, ``i_1``, ``i_2`` etc. Unlike
Clifford algebras, these numbers commute, i.e. ``i_1i_2=i_2i_1``.

## Algebra of the multicomplex numbers

### Recursive Definition

Multicomplex numbers are built recursively from real numbers:

```math
\mathbb{C}_0 = \mathbb{R}
```
```math
\mathbb{C}_n = \mathbb{C}_{n-1} + i_n \mathbb{C}_{n-1}
```

This means:
- ``\mathbb{C}_0``: Real numbers (1 component)
- ``\mathbb{C}_1 = \mathbb{R} + i_1\mathbb{R}``: Complex numbers (2 components)
- ``\mathbb{C}_2 = \mathbb{C}_1 + i_2\mathbb{C}_1``: Bicomplex numbers (4 components)
- ``\mathbb{C}_3 = \mathbb{C}_2 + i_3\mathbb{C}_2``: Tricomplex numbers (8 components)
- ``\mathbb{C}_n``: Has ``2^n`` real components

### Fundamental Properties

Each imaginary unit squares to negative one:
```math
i_k^2 = -1 \quad \text{for all } k \geq 1
```

In code:
```@repl background
using MulticomplexNumbers

im1 * im1
im2 * im2
```

Crucially, imaginary units **commute** with each other:
```math
i_j i_k = i_k i_j \quad \text{for all } j, k
```

This commutativity distinguishes multicomplex numbers from quaternions and other hypercomplex systems where imaginary units anticommute.

```@repl background
# Commutativity
im1 * im2 == im2 * im1
```

Interestingly, products of *distinct* imaginary units square to **positive one**:
```math
(i_j i_k)^2 = i_j i_k i_j i_k = i_j^2 i_k^2 = (-1)(-1) = +1 \quad \text{(when } j \neq k \text{)}
```

```@repl background
# Hyperbolic units
(im1 * im2)^2
```

This reveals an important property: products like ``i_1 i_2`` are **hyperbolic units** (they square to +1), not elliptic units (which square to -1). This is why ``\mathbb{C}_n`` for ``n \geq 2`` contains zero divisors.

### Component Representation

A general multicomplex number ``m \in \mathbb{C}_n`` can be written as:

For ``n = 1`` (complex):
```math
m = a + b \cdot i_1
```

For ``n = 2`` (bicomplex):
```math
m = a + b \cdot i_1 + c \cdot i_2 + d \cdot i_1 i_2
```

For ``n = 3`` (tricomplex):
```math
m = a + b \cdot i_1 + c \cdot i_2 + d \cdot i_1 i_2 + e \cdot i_3 + f \cdot i_1 i_3 + g \cdot i_2 i_3 + h \cdot i_1 i_2 i_3
```

In general, ``\mathbb{C}_n`` has ``2^n`` basis elements formed by all possible products of the imaginary units ``\{1, i_1, i_2, \ldots, i_n\}`` and their combinations.

### Recursive Structure

Any multicomplex number ``m \in \mathbb{C}_n`` can be uniquely written as:
```math
m = a + i_n \cdot b
```
where ``a, b \in \mathbb{C}_{n-1}`` are called the "real" and "imaginary" parts (with respect to ``i_n``).

This recursive structure is central to the implementation: the `real` and `imag` functions return multicomplex numbers of order ``n-1``.

```@repl background
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2

# Recursive structure: z = real(z) + im2 * imag(z)
real(z)   # Order 1 (with respect to im2)
imag(z)   # Order 1
```

## Matrix Representations

### The Matrix Homomorphism

Every multicomplex number can be faithfully represented as a real matrix. This representation is a ring homomorphism that preserves addition and multiplication.

For ``\mathbb{C}_1`` (complex numbers), the familiar 2×2 representation:
```math
a + b \cdot i_1 \mapsto \begin{pmatrix} a & -b \\ b & a \end{pmatrix}
```

```@repl background
using MulticomplexNumbers

z = 3.0 + 4.0*im1
matrep(z)
```

For ``\mathbb{C}_2`` (bicomplex numbers), a 4×4 matrix where the structure is built recursively:
```math
(a + b \cdot i_1) + (c + d \cdot i_1) \cdot i_2 \mapsto \begin{pmatrix} M(a + b \cdot i_1) & -M(c + d \cdot i_1) \\ M(c + d \cdot i_1) & M(a + b \cdot i_1) \end{pmatrix}
```

where ``M(\cdot)`` denotes the matrix representation.

```@repl background
w = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2
matrep(w)
```

### Recursive Construction

For ``m = a + i_n \cdot b \in \mathbb{C}_n`` with ``a, b \in \mathbb{C}_{n-1}``:
```math
M(m) = \begin{pmatrix} M(a) & -M(b) \\ M(b) & M(a) \end{pmatrix}
```

This ``2^n \times 2^n`` matrix has the property that the first column contains all ``2^n`` components of the multicomplex number.

### Why Matrix Representation?

The matrix representation enables:
1. **Multiplication**: ``m_1 \times m_2`` computed via ``M(m_1) \times M(m_2)``
2. **Division**: ``m_1 / m_2`` computed via ``M(m_1) \times M(m_2)^{-1}``
3. **Powers**: ``m^p`` computed via matrix power
4. **Transcendental functions**: ``\exp(m)``, ``\log(m)``, ``\sqrt{m}`` via matrix functions

The result is extracted from the first column of the resulting matrix.

```@repl background
using MulticomplexNumbers
using LinearAlgebra

z = 1.0 + 0.5*im1

# Verify: multiplication via matrices
w = 2.0 + 0.3*im1
z * w

# Compare to matrix multiplication
M_z = matrep(z)
M_w = matrep(w)
M_product = M_z * M_w
M_product[:, 1]  # First column contains the components

# Transcendental functions use matrix functions
exp(z)
exp(M_z)[:, 1]  # Same result
```

### Zero Divisors

Unlike complex numbers, bicomplex and higher multicomplex numbers have **zero divisors**: non-zero elements whose product is zero.

For example, in ``\mathbb{C}_2``:
```math
(1 + i_1 i_2)(1 - i_1 i_2) = 1 - (i_1 i_2)^2 = 1 - 1 = 0
```

This means:
- Division by some non-zero elements is undefined (singular matrix)
- ``\mathbb{C}_n`` for ``n \geq 2`` is not a division algebra

## Arithmetic

### Addition and Subtraction

Component-wise operations:
```math
(a + b \cdot i_1) + (c + d \cdot i_1) = (a + c) + (b + d) \cdot i_1
```

This extends naturally to all orders.

```@repl background
using MulticomplexNumbers

a = 1.0 + 2.0*im1
b = 3.0 + 4.0*im1
a + b
```

### Multiplication

MulticomplexNumbers.jl uses a **recursive algorithm** for multiplication that is more efficient than the matrix representation approach.

For ``w = u + v \cdot i_n`` and ``z = x + y \cdot i_n`` where ``u,v,x,y \in \mathbb{C}_{n-1}``:
```math
w \cdot z = (u \cdot x - v \cdot y) + (u \cdot y + v \cdot x) \cdot i_n
```

This recursively reduces to real number multiplication at the base case.

```@repl background
a = 1.0 + 2.0*im1
b = 3.0 + 4.0*im1
a * b
```

The recursive approach is more efficient than matrix multiplication for most cases, especially for lower orders.

**References:**
- Implementation follows the recursive definition from NIST report: https://doi.org/10.6028/jres.126.033

### Division

Division is computed using a **recursive conjugate-folding algorithm** that avoids matrix operations.

The **fold operation** is defined as:
```math
\text{fold}(m) = m \cdot \overline{m} \in \mathbb{C}_{n-1}
```

This reduces the order by one. For division ``a / b``:
```math
\frac{a}{b} = \frac{a \cdot \overline{b}}{\text{fold}(b)}
```

Each fold reduces the denominator's order by 1. After ``n`` folds, the denominator becomes a real number, enabling simple scalar division.

For **abient numbers** (those whose repeated folding produces zero), the algorithm falls back to matrix inversion.

**References:**
- Recursive division algorithm from NIST report: https://doi.org/10.6028/jres.126.033

### Conjugation

Multicomplex conjugation inverts the sign of the highest imaginary component:
```math
\overline{m} = \overline{a + i_n \cdot b} = a - i_n \cdot b
```

where ``a, b \in \mathbb{C}_{n-1}``.

Note: This is different from complex conjugation applied to all imaginary units.

```@repl background
using MulticomplexNumbers

z = 1.0 + 2.0*im1 + 3.0*im2 + 4.0*im1*im2
conj(z)  # Negates the im2 component
```

### Norm

The Euclidean norm is defined as:
```math
\|m\| = \sqrt{\sum_{k=1}^{2^n} c_k^2}
```

where ``c_k`` are the real components of the multicomplex number.

```@repl background
using MulticomplexNumbers
using LinearAlgebra

z = 3.0 + 4.0*im1
norm(z)  # sqrt(3² + 4²) = 5
abs(z)   # Same as norm
```

### Trigonometric and Hyperbolic Functions

MulticomplexNumbers.jl implements trigonometric and hyperbolic functions using **recursive algorithms** that are more efficient than matrix-based methods.

For a multicomplex number ``x = x_0 + i_n x_1`` where ``x_0, x_1 \in \mathbb{C}_{n-1}``:

**Trigonometric functions:**
```math
\sin(x) = \sin(x_0)\cosh(x_1) + i_n \cos(x_0)\sinh(x_1)
```
```math
\cos(x) = \cos(x_0)\cosh(x_1) - i_n \sin(x_0)\sinh(x_1)
```

**Hyperbolic functions:**
```math
\sinh(x) = \sinh(x_0)\cos(x_1) + i_n \cosh(x_0)\sin(x_1)
```
```math
\cosh(x) = \cosh(x_0)\cos(x_1) + i_n \sinh(x_0)\sin(x_1)
```

These formulas recursively reduce to real number operations at the base case (``n=0``). For order 1 multicomplex numbers (standard complex numbers), the implementation uses Julia's optimized complex arithmetic directly.

Other functions are derived from these:
- ``\tan(x) = \sin(x) / \cos(x)``
- ``\tanh(x) = \sinh(x) / \cosh(x)``
- ``\sec(x) = 1 / \cos(x)``, ``\csc(x) = 1 / \sin(x)``
- ``\sech(x) = 1 / \cosh(x)``, ``\csch(x) = 1 / \sinh(x)``

Inverse trigonometric and hyperbolic functions currently use matrix representations, as recursive formulas are not provided in the reference literature.

```@repl background
using MulticomplexNumbers

z = 0.5 + 0.3*im1
sin(z)
cos(z)^2 + sin(z)^2  # Pythagorean identity holds

w = 1.0 + 0.5*im1 + 0.3*im2
sinh(w)
cosh(w)^2 - sinh(w)^2  # Hyperbolic identity holds
```

**References:**
- Recursive formulas from NIST report: https://doi.org/10.6028/jres.126.033 (page 15)

## Applications

### Numerical Differentiation

The primary application of multicomplex numbers is **high-order numerical differentiation**. The multicomplex step method extends the complex step derivative to compute higher-order derivatives with machine precision accuracy.

For a function ``f`` extended to multicomplex arguments:
```math
f(x + h \cdot i_1) = f(x) + h \cdot f'(x) \cdot i_1 + O(h^2)
```

The first derivative is extracted from the ``i_1`` component:
```math
f'(x) \approx \frac{\text{Im}_1[f(x + h \cdot i_1)]}{h}
```

This avoids subtractive cancellation errors that plague finite difference methods.

For second derivatives using bicomplex numbers:
```math
f''(x) \approx \frac{\text{component}_{i_1 i_2}[f(x + h \cdot i_1 + h \cdot i_2)]}{h^2}
```

### Signal Processing

Multicomplex Fourier transforms enable processing of signals with multiple independent phase components, useful in certain NMR spectroscopy applications.

## References

1. NIST report on multicomplex algebra: [https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf](https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf)
2. Casado & Hewson (2020): Algorithm 1008: Multicomplex Number Class for Matlab. ACM Trans Math Softw. [http://dx.doi.org/10.1145/3378542](http://dx.doi.org/10.1145/3378542)
3. NIST C++ implementation: [https://github.com/usnistgov/multicomplex](https://github.com/usnistgov/multicomplex)