"invert if y<0"
Base.flipsign(x::Multicomplex, y::Real) = ifelse(signbit(y), -x, x)


###############
# comparisons #
###############
"equality if all components are equal"
Base.:(==)(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = (a.value == b.value)
Base.isequal(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = isequal(a.value, b.value)

"compute hash from array of components"
Base.hash(m::Multicomplex, h::UInt) = hash(m.value, h)


###############
# conjugation #
###############
@doc raw"""
    conj(m::Multicomplex{T,N})

Multicomplex conjugation in ℂₙ is defined by inverting the sign of the
imaginary component: ``conj(mₙ = aₙ₋₁ + iₙbₙ₋₁) = mₙ = aₙ₋₁ - iₙbₙ₋₁``.
"""
Base.conj(m::Multicomplex) = Multicomplex(real(m), -imag(m))


###################
# absolute values #
###################
Base.abs2(m::Multicomplex{T,N}) where {T,N} = sum(x -> x^2, m.value)
Base.abs(m::Multicomplex{T,N}) where {T,N} = sqrt(abs2(m))


####################
# basic arithmetic #
####################
"unary plus"
Base.:(+)(a::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(a.value)

"unary minus"
Base.:(-)(a::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(-a.value)

"addition of multicomplex numbers: calculated as sum of components"
Base.:(+)(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(a.value + b.value)
Base.:(+)(a::Multicomplex{T,N}, b::Real) where {T,N} = a + Multicomplex(b)
Base.:(+)(a::Real, b::Multicomplex{T,N}) where {T,N} = Multicomplex(a) + b

"subtraction of multicomplex numbers: calculated as difference of components"
Base.:(-)(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(a.value - b.value)
Base.:(-)(a::Multicomplex{T,N}, b::Real) where {T,N} = a - Multicomplex(b)
Base.:(-)(a::Real, b::Multicomplex{T,N}) where {T,N} = Multicomplex(a) - b

"multiplication by a scalar"
Base.:(*)(a::Real, m::Multicomplex{T,N,C}) where {T,N,C} = Multicomplex{N}(a * m.value)
Base.:(*)(m::Multicomplex{T,N,C}, a::Real) where {T,N,C} = Multicomplex{N}(a * m.value)

"division by a scalar"
Base.:(/)(m::Multicomplex{T,N,C}, a::Real) where {T,N,C} = Multicomplex{N}(m.value / a)

##################
# multiplication #
##################
@doc raw"""
    *(a, b)

Multicomplex multiplication using recursive definition.

For two multicomplex numbers ``w,z \in \mathbb{C}_{n+1}``, written as
``w = u + v \cdot i_{n+1}`` and ``z = x + y \cdot i_{n+1}`` where ``u,v,x,y \in \mathbb{C}_n``,
the multiplication rule is:

``w \cdot z = (u \cdot x - v \cdot y) + (u \cdot y + v \cdot x) \cdot i_{n+1}``

This is more efficient than the matrix representation approach for most cases.
"""
# Base case: N=0 (real numbers wrapped in Multicomplex)
function Base.:(*)(a::Multicomplex{T,0,1}, b::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(a.value[1] * b.value[1]))
end

# Recursive case: N≥1
# w = u + v*i_N, z = x + y*i_N
# w*z = (u*x - v*y) + (u*y + v*x)*i_N
function Base.:(*)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    # Extract real and imaginary parts (each of order N-1)
    u = real(a)  # Multicomplex{N-1}
    v = imag(a)  # Multicomplex{N-1}
    x = real(b)  # Multicomplex{N-1}
    y = imag(b)  # Multicomplex{N-1}

    # Compute products recursively
    ux = u * x
    vy = v * y
    uy = u * y
    vx = v * x

    # Combine: (ux - vy) + (uy + vx)*i_N
    real_part = ux - vy
    imag_part = uy + vx

    Multicomplex(real_part, imag_part)
end



##################
# powers         #
##################
"""
    ^(a,b::Real)
Multicomplex multiplication via the matrix representation.
"""
# Optimized case: N=1 (standard complex arithmetic)
function Base.:(^)(a::Multicomplex{T,1,2}, b::Real) where {T}
    Multicomplex(Complex(a.value[1], a.value[2])^b)
end

# General case: use matrix representation
function Base.:(^)(a::Multicomplex{T,N,C}, b::Real) where {T,N,C}
    Multicomplex{N}(SMatrix{C,C}(real(matrep(a)^b))[SVector{C}(SOneTo(C))])
end

# Optimized case: N=1 (standard complex arithmetic) - integer specialization
function Base.:(^)(a::Multicomplex{T,1,2}, b::Integer) where {T}
    Multicomplex(Complex(a.value[1], a.value[2])^b)
end

# General case: integer specialization needed to avoid method ambiguity
function Base.:(^)(a::Multicomplex{T,N,C}, b::Integer) where {T,N,C}
    Multicomplex{N}(SMatrix{C,C}(real(matrep(a)^b))[SVector{C}(SOneTo(C))])
end



############
# division #
############
@doc raw"""
    inv(m::Multicomplex)

Compute the multiplicative inverse of a multicomplex number using a recursive
conjugate-folding algorithm.

For a non-abient multicomplex number ``m \in \mathbb{C}_n``:
```math
m^{-1} = \overline{m} \cdot (\text{fold}(m))^{-1}
```

where ``\overline{m}`` is the conjugate and ``\text{fold}(m) = m \cdot \overline{m} \in \mathbb{C}_{n-1}``.

This recursively reduces the order until we reach a real number, where inversion
is simply taking the reciprocal.

If `m` is abient (i.e., folds to zero), falls back to matrix representation
which will throw `LinearAlgebra.SingularException`.
"""
# Base case: N=0 (real numbers)
function Base.inv(m::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(inv(m.value[1])))
end

# Recursive case: N≥1
# inv(m) = conj(m) * inv(fold(m))
function Base.inv(m::Multicomplex{T,N,C}) where {T,N,C}
    folded = fold(m)
    # Check if fold result is approximately zero (abient case)
    # If so, fall back to matrix representation
    if iszero(folded)
        # Use matrix representation for abient numbers
        return Multicomplex{N}((I(C)/matrep(m))[SVector{C}(SOneTo(C))])
    end
    # Recursive inversion: inv(m) = conj(m) * inv(fold(m))
    conj(m) * inv(folded)
end

@doc raw"""
    /(a, b)

Multicomplex division using recursive conjugate-folding algorithm.

The algorithm works by repeatedly multiplying numerator and denominator by the
conjugate of the denominator:
```math
\frac{a}{b} = \frac{a \cdot \overline{b}}{b \cdot \overline{b}} = \frac{a \cdot \overline{b}}{\text{fold}(b)}
```

Each fold operation reduces the order of the denominator by 1. After ``n`` folds,
the denominator becomes a real number, and division is simply scalar division.

For abient divisors (where fold eventually produces zero), falls back to matrix
representation which will throw `LinearAlgebra.SingularException`.
"""
function Base.:(/)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    a * inv(b)
end



################
# exponentials #
################
@doc raw"""
    exp(m::Multicomplex)

Multicomplex exponential using recursive algorithm.

For a multicomplex number ``a = x_0 + i_n x_1`` where ``x_0, x_1 \in \mathbb{C}_{n-1}``:

```math
\exp(a) = \exp(x_0) \cos(x_1) + i_n \exp(x_0) \sin(x_1)
```

Equivalently: ``\exp(x_0 + i_n x_1) = \exp(x_0) [\cos(x_1) + i_n \sin(x_1)]``

This recursive approach is more efficient than matrix representation.

# References
- NIST report equations (44-45): https://doi.org/10.6028/jres.126.033
"""
# Base case: real numbers
function Base.exp(m::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(exp(m.value[1])))
end

# Optimized case: N=1 (standard complex arithmetic)
function Base.exp(m::Multicomplex{T,1,2}) where {T}
    Multicomplex(exp(Complex(m.value[1], m.value[2])))
end

# Recursive case: N≥2
# Always use equation 44 - simpler and more reliable
function Base.exp(m::Multicomplex{T,N,C}) where {T,N,C}
    # Equation 44: exp(x₀ + i_n x₁) = exp(x₀)[cos(x₁) + i_n sin(x₁)]
    x0 = real(m)  # Multicomplex{N-1}
    x1 = imag(m)  # Multicomplex{N-1}

    exp_x0 = exp(x0)
    y0 = exp_x0 * cos(x1)
    y1 = exp_x0 * sin(x1)

    Multicomplex(y0, y1)
end


################
# logarithm    #
################
"""
    log(m)
Multicomplex logarithm via the matrix representation.
"""
# Optimized case: N=1 (standard complex arithmetic)
function Base.log(m::Multicomplex{T,1,2}) where {T}
    Multicomplex(log(Complex(m.value[1], m.value[2])))
end

# General case: use matrix representation
function Base.log(m::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}(SMatrix{C,C}(log(matrep(m)))[SVector{C}(SOneTo(C))])
end


################
# square root  #
################
"""
    sqrt(m)
Multicomplex square root via the matrix representation.
"""
# Optimized case: N=1 (standard complex arithmetic)
function Base.sqrt(m::Multicomplex{T,1,2}) where {T}
    Multicomplex(sqrt(Complex(m.value[1], m.value[2])))
end

# General case: use matrix representation
function Base.sqrt(m::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}(SMatrix{C,C}(sqrt(matrep(m)))[SVector{C}(SOneTo(C))])
end


################
# fold operator#
################
@doc raw"""
    fold(m::Multicomplex{T,N})

Folding operation that multiplies a multicomplex number by its conjugate.

With w ∈ ℂₙ₊₁, fold(w) = w × conj(w) = conj(w) × w ∈ ℂₙ.

The fold operation reduces the order by 1:
- fold(ℂₙ) → ℂₙ₋₁
- fold(ℂ₁) → ℂ₀ (real numbers)
- fold(ℂ₀) → ℂ₀ (squaring)

# Examples
```julia
julia> fold(1 + im1)
Multicomplex(2.0)

julia> fold(1 + im1*im2)  # (1 + i₁i₂)(1 - i₁i₂) = 0
Multicomplex(0.0, 0.0)
```

# References
Sometimes the fold of a nonzero number turns out to be zero; for example,
(1 + i₁i₂)(1 - i₁i₂) = 0. This illustrates that the norm is not always
preserved under multiplication, thus the multicomplex numbers are not a
composition algebra.
"""
function fold(m::Multicomplex{T,0,1}) where {T}
    # For real numbers (ℂ₀), fold is just squaring
    Multicomplex{0}(SVector(m.value[1]^2))
end

function fold(m::Multicomplex{T,1,2}) where {T}
    # For complex numbers (ℂ₁), fold gives |z|² as a real number (ℂ₀)
    result = m * conj(m)
    Multicomplex{0}(SVector(real(result)))
end

function fold(m::Multicomplex{T,N,C}) where {T,N,C}
    # For higher orders (ℂₙ where n≥2), fold reduces order by 1
    result = m * conj(m)
    real(result)  # Returns Multicomplex{N-1}
end

@doc raw"""
    isabient(m::Multicomplex; atol=0, rtol=atol>0 ? 0 : √eps)

Determine if a multicomplex number is abient.

A multicomplex number is abient (from a Latin participle meaning "going away")
if after sufficiently many foldings it becomes zero. After N foldings, an ℂₙ
multicomplex number becomes a real number; if it turns out to be zero, the
multicomplex number is abient.

# Arguments
- `m::Multicomplex`: The multicomplex number to test
- `atol::Real=0`: Absolute tolerance for zero comparison
- `rtol::Real`: Relative tolerance (defaults to √eps if atol=0)

# Examples
```julia
julia> isabient(1 + im1)
false

julia> isabient(1 + im1*im2)  # (1 + i₁i₂) is abient
true
```

# References
Not to be confused with "ambient". This concept arises because multicomplex
numbers are not a composition algebra - the norm is not always preserved
under multiplication.
"""
function isabient(m::Multicomplex{T,N,C}; atol::Real=0, rtol::Real=Base.rtoldefault(T)) where {T,N,C}
    # Fold N times to reduce to a real number (ℂ₀)
    result = m
    for _ in 1:N
        result = fold(result)
    end
    # After N folds, result is Multicomplex{0} - check if it's approximately zero
    # For Multicomplex{0}, we check the single value component
    return isapprox(result.value[1], zero(T); atol=atol, rtol=rtol)
end


#########################
# Trigonometric functions
#########################

@doc raw"""
    sin(m::Multicomplex)

Multicomplex sine function using recursive algorithm.

For a multicomplex number ``x = x_0 + i_n x_1`` where ``x_0, x_1 \in \mathbb{C}_{n-1}``:

```math
\sin(x) = \sin(x_0)\cosh(x_1) + i_n \cos(x_0)\sinh(x_1)
```

This recursive approach is more efficient than matrix representation.

# References
- NIST report on multicomplex algebra: https://doi.org/10.6028/jres.126.033
"""
# Base case: real numbers
function Base.sin(m::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(sin(m.value[1])))
end

# Optimized case: N=1 (standard complex arithmetic)
function Base.sin(m::Multicomplex{T,1,2}) where {T}
    z = Complex(m.value[1], m.value[2])
    Multicomplex(sin(z))
end

# Recursive case: N≥2
function Base.sin(m::Multicomplex{T,N,C}) where {T,N,C}
    x0 = real(m)  # Multicomplex{N-1}
    x1 = imag(m)  # Multicomplex{N-1}

    # sin(x) = sin(x0)*cosh(x1) + i_N*cos(x0)*sinh(x1)
    y0 = sin(x0) * cosh(x1)
    y1 = cos(x0) * sinh(x1)

    Multicomplex(y0, y1)
end

@doc raw"""
    cos(m::Multicomplex)

Multicomplex cosine function using recursive algorithm.

For a multicomplex number ``x = x_0 + i_n x_1`` where ``x_0, x_1 \in \mathbb{C}_{n-1}``:

```math
\cos(x) = \cos(x_0)\cosh(x_1) - i_n \sin(x_0)\sinh(x_1)
```

This recursive approach is more efficient than matrix representation.

# References
- NIST report on multicomplex algebra: https://doi.org/10.6028/jres.126.033
"""
# Base case: real numbers
function Base.cos(m::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(cos(m.value[1])))
end

# Optimized case: N=1 (standard complex arithmetic)
function Base.cos(m::Multicomplex{T,1,2}) where {T}
    z = Complex(m.value[1], m.value[2])
    Multicomplex(cos(z))
end

# Recursive case: N≥2
function Base.cos(m::Multicomplex{T,N,C}) where {T,N,C}
    x0 = real(m)  # Multicomplex{N-1}
    x1 = imag(m)  # Multicomplex{N-1}

    # cos(x) = cos(x0)*cosh(x1) - i_N*sin(x0)*sinh(x1)
    y0 = cos(x0) * cosh(x1)
    y1 = -sin(x0) * sinh(x1)

    Multicomplex(y0, y1)
end

"""
    tan(m::Multicomplex)

Multicomplex tangent function computed as sin/cos.
"""
Base.tan(m::Multicomplex) = sin(m) / cos(m)

"""
    cot(m::Multicomplex)

Multicomplex cotangent function computed as cos/sin.
"""
Base.cot(m::Multicomplex) = cos(m) / sin(m)

"""
    sec(m::Multicomplex)

Multicomplex secant function computed as 1/cos.
"""
Base.sec(m::Multicomplex) = inv(cos(m))

"""
    csc(m::Multicomplex)

Multicomplex cosecant function computed as 1/sin.
"""
Base.csc(m::Multicomplex) = inv(sin(m))


############################
# Hyperbolic functions
############################

@doc raw"""
    sinh(m::Multicomplex)

Multicomplex hyperbolic sine function using recursive algorithm.

For a multicomplex number ``x = x_0 + i_n x_1`` where ``x_0, x_1 \in \mathbb{C}_{n-1}``:

```math
\sinh(x) = \sinh(x_0)\cos(x_1) + i_n \cosh(x_0)\sin(x_1)
```

This recursive approach is more efficient than matrix representation.

# References
- NIST report on multicomplex algebra: https://doi.org/10.6028/jres.126.033
"""
# Base case: real numbers
function Base.sinh(m::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(sinh(m.value[1])))
end

# Optimized case: N=1 (standard complex arithmetic)
function Base.sinh(m::Multicomplex{T,1,2}) where {T}
    z = Complex(m.value[1], m.value[2])
    Multicomplex(sinh(z))
end

# Recursive case: N≥2
function Base.sinh(m::Multicomplex{T,N,C}) where {T,N,C}
    x0 = real(m)  # Multicomplex{N-1}
    x1 = imag(m)  # Multicomplex{N-1}

    # sinh(x) = sinh(x0)*cos(x1) + i_N*cosh(x0)*sin(x1)
    y0 = sinh(x0) * cos(x1)
    y1 = cosh(x0) * sin(x1)

    Multicomplex(y0, y1)
end

@doc raw"""
    cosh(m::Multicomplex)

Multicomplex hyperbolic cosine function using recursive algorithm.

For a multicomplex number ``x = x_0 + i_n x_1`` where ``x_0, x_1 \in \mathbb{C}_{n-1}``:

```math
\cosh(x) = \cosh(x_0)\cos(x_1) + i_n \sinh(x_0)\sin(x_1)
```

This recursive approach is more efficient than matrix representation.

# References
- NIST report on multicomplex algebra: https://doi.org/10.6028/jres.126.033
"""
# Base case: real numbers
function Base.cosh(m::Multicomplex{T,0,1}) where {T}
    Multicomplex{0}(SVector(cosh(m.value[1])))
end

# Optimized case: N=1 (standard complex arithmetic)
function Base.cosh(m::Multicomplex{T,1,2}) where {T}
    z = Complex(m.value[1], m.value[2])
    Multicomplex(cosh(z))
end

# Recursive case: N≥2
function Base.cosh(m::Multicomplex{T,N,C}) where {T,N,C}
    x0 = real(m)  # Multicomplex{N-1}
    x1 = imag(m)  # Multicomplex{N-1}

    # cosh(x) = cosh(x0)*cos(x1) + i_N*sinh(x0)*sin(x1)
    y0 = cosh(x0) * cos(x1)
    y1 = sinh(x0) * sin(x1)

    Multicomplex(y0, y1)
end

"""
    tanh(m::Multicomplex)

Multicomplex hyperbolic tangent function computed as sinh/cosh.
"""
Base.tanh(m::Multicomplex) = sinh(m) / cosh(m)

"""
    coth(m::Multicomplex)

Multicomplex hyperbolic cotangent function computed as cosh/sinh.
"""
Base.coth(m::Multicomplex) = cosh(m) / sinh(m)

"""
    sech(m::Multicomplex)

Multicomplex hyperbolic secant function computed as 1/cosh.
"""
Base.sech(m::Multicomplex) = inv(cosh(m))

"""
    csch(m::Multicomplex)

Multicomplex hyperbolic cosecant function computed as 1/sinh.
"""
Base.csch(m::Multicomplex) = inv(sinh(m))


####################################
# Inverse trigonometric functions
####################################

# Helper to get matrix in the right format for matrix-based trig functions
# Julia 1.11 requires mutable matrices for some operations, Julia 1.12+ works with immutable
@inline _matfunc_matrix(mat) = @static if VERSION >= v"1.12"
    mat  # Use immutable SMatrix directly in Julia 1.12+
else
    Matrix(mat)  # Convert to mutable Matrix for Julia 1.11
end

"""
    asin(m::Multicomplex)

Multicomplex arcsine function via the matrix representation.
"""
function Base.asin(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(asin(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    acos(m::Multicomplex)

Multicomplex arccosine function via the matrix representation.
"""
function Base.acos(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(acos(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    atan(m::Multicomplex)

Multicomplex arctangent function via the matrix representation.
"""
function Base.atan(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(atan(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    acot(m::Multicomplex)

Multicomplex arccotangent function via the matrix representation.
"""
function Base.acot(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(acot(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    asec(m::Multicomplex)

Multicomplex arcsecant function via the matrix representation.
"""
function Base.asec(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(asec(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    acsc(m::Multicomplex)

Multicomplex arccosecant function via the matrix representation.
"""
function Base.acsc(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(acsc(_matfunc_matrix(matrep(m)))[:, 1])))
end


####################################
# Inverse hyperbolic functions
####################################

"""
    asinh(m::Multicomplex)

Multicomplex inverse hyperbolic sine function via the matrix representation.
"""
function Base.asinh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(asinh(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    acosh(m::Multicomplex)

Multicomplex inverse hyperbolic cosine function via the matrix representation.
"""
function Base.acosh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(acosh(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    atanh(m::Multicomplex)

Multicomplex inverse hyperbolic tangent function via the matrix representation.
"""
function Base.atanh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(atanh(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    acoth(m::Multicomplex)

Multicomplex inverse hyperbolic cotangent function via the matrix representation.
"""
function Base.acoth(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(acoth(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    asech(m::Multicomplex)

Multicomplex inverse hyperbolic secant function via the matrix representation.
"""
function Base.asech(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(asech(_matfunc_matrix(matrep(m)))[:, 1])))
end

"""
    acsch(m::Multicomplex)

Multicomplex inverse hyperbolic cosecant function via the matrix representation.
"""
function Base.acsch(m::Multicomplex{T,N,C}) where {T,N,C}
    # Inverse functions can return complex values, so extract real parts
    Multicomplex{N}(SVector{C}(real.(acsch(_matfunc_matrix(matrep(m)))[:, 1])))
end