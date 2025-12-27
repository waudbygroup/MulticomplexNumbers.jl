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
Base.conj(m::Multicomplex) = Multicomplex(real(m),-imag(m))


###################
# absolute values #
###################
Base.abs2(m::Multicomplex{T,N}) where {T,N} = sum(x->x^2, m.value)
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
Base.:(*)(a::Real, m::Multicomplex{T,N,C}) where {T,N,C} = Multicomplex{N}(a*m.value)
Base.:(*)(m::Multicomplex{T,N,C}, a::Real) where {T,N,C} = Multicomplex{N}(a*m.value)

"division by a scalar"
Base.:(/)(m::Multicomplex{T,N,C}, a::Real) where {T,N,C} = Multicomplex{N}(m.value / a)

##################
# multiplication #
##################
"""
    *(a,b)
Multicomplex multiplication via the matrix representation.
"""
function Base.:(*)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}((matrep(a) * matrep(b))[SVector{C}(SOneTo(C))])
end



##################
# powers         #
##################
"""
    ^(a,b::Real)
Multicomplex multiplication via the matrix representation.
"""
function Base.:(^)(a::Multicomplex{T,N,C}, b::Real) where {T,N,C}
    Multicomplex{N}(SMatrix{C,C}(real(matrep(a) ^ b))[SVector{C}(SOneTo(C))])
end
# integer specialisation needed to avoid method ambiguity
function Base.:(^)(a::Multicomplex{T,N,C}, b::Integer) where {T,N,C}
    Multicomplex{N}(SMatrix{C,C}(real(matrep(a) ^ b))[SVector{C}(SOneTo(C))])
end



############
# division #
############
"""
    /(a,b)
Multicomplex division via the matrix representation.

Throws LinearAlgebra.SingularException if zero divisors exist (det b = 0).
"""
function Base.:(/)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}((matrep(a) / matrep(b))[SVector{C}(SOneTo(C))])
end

function Base.inv(m::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}((I / matrep(m))[SVector{C}(SOneTo(C))])
end
Base.inv(m::Multicomplex{<:Integer}) = inv(float(m))



################
# exponentials #
################
"""
    exp(m)
Multicomplex exponential via the matrix representation.
"""
function Base.exp(m::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}(exp(matrep(m))[SVector{C}(SOneTo(C))])
end


################
# logarithm    #
################
"""
    log(m)
Multicomplex logarithm via the matrix representation.
"""
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

The fold of a ℂ₀ number (real number) is simply its square.

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
function fold(m::Multicomplex{T,N,C}) where {T,N,C}
    m * conj(m)
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
    # Fold N times to reduce to a real number
    result = m
    for _ in 1:N
        result = fold(result)
    end
    # Check if the resulting real number is approximately zero
    return isapprox(realest(result), zero(T); atol=atol, rtol=rtol)
end


#########################
# Trigonometric functions
#########################

"""
    sin(m::Multicomplex)

Multicomplex sine function via the matrix representation.
"""
function Base.sin(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(sin(mat))[SVector{C}(SOneTo(C))])
end

"""
    cos(m::Multicomplex)

Multicomplex cosine function via the matrix representation.
"""
function Base.cos(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(cos(mat))[SVector{C}(SOneTo(C))])
end

"""
    tan(m::Multicomplex)

Multicomplex tangent function via the matrix representation.
"""
function Base.tan(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(tan(mat))[SVector{C}(SOneTo(C))])
end

"""
    cot(m::Multicomplex)

Multicomplex cotangent function via the matrix representation.
"""
function Base.cot(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(cot(mat))[SVector{C}(SOneTo(C))])
end

"""
    sec(m::Multicomplex)

Multicomplex secant function via the matrix representation.
"""
function Base.sec(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(sec(mat))[SVector{C}(SOneTo(C))])
end

"""
    csc(m::Multicomplex)

Multicomplex cosecant function via the matrix representation.
"""
function Base.csc(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(csc(mat))[SVector{C}(SOneTo(C))])
end


############################
# Hyperbolic functions
############################

"""
    sinh(m::Multicomplex)

Multicomplex hyperbolic sine function via the matrix representation.
"""
function Base.sinh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(sinh(mat))[SVector{C}(SOneTo(C))])
end

"""
    cosh(m::Multicomplex)

Multicomplex hyperbolic cosine function via the matrix representation.
"""
function Base.cosh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(cosh(mat))[SVector{C}(SOneTo(C))])
end

"""
    tanh(m::Multicomplex)

Multicomplex hyperbolic tangent function via the matrix representation.
"""
function Base.tanh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(tanh(mat))[SVector{C}(SOneTo(C))])
end

"""
    coth(m::Multicomplex)

Multicomplex hyperbolic cotangent function via the matrix representation.
"""
function Base.coth(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(coth(mat))[SVector{C}(SOneTo(C))])
end

"""
    sech(m::Multicomplex)

Multicomplex hyperbolic secant function via the matrix representation.
"""
function Base.sech(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(sech(mat))[SVector{C}(SOneTo(C))])
end

"""
    csch(m::Multicomplex)

Multicomplex hyperbolic cosecant function via the matrix representation.
"""
function Base.csch(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    Multicomplex{N}(real.(csch(mat))[SVector{C}(SOneTo(C))])
end


####################################
# Inverse trigonometric functions
####################################

"""
    asin(m::Multicomplex)

Multicomplex arcsine function via the matrix representation.
"""
function Base.asin(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = asin(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    acos(m::Multicomplex)

Multicomplex arccosine function via the matrix representation.
"""
function Base.acos(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = acos(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    atan(m::Multicomplex)

Multicomplex arctangent function via the matrix representation.
"""
function Base.atan(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = atan(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    acot(m::Multicomplex)

Multicomplex arccotangent function via the matrix representation.
"""
function Base.acot(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = acot(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    asec(m::Multicomplex)

Multicomplex arcsecant function via the matrix representation.
"""
function Base.asec(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = asec(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    acsc(m::Multicomplex)

Multicomplex arccosecant function via the matrix representation.
"""
function Base.acsc(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = acsc(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end


####################################
# Inverse hyperbolic functions
####################################

"""
    asinh(m::Multicomplex)

Multicomplex inverse hyperbolic sine function via the matrix representation.
"""
function Base.asinh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = asinh(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    acosh(m::Multicomplex)

Multicomplex inverse hyperbolic cosine function via the matrix representation.
"""
function Base.acosh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = acosh(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    atanh(m::Multicomplex)

Multicomplex inverse hyperbolic tangent function via the matrix representation.
"""
function Base.atanh(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = atanh(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    acoth(m::Multicomplex)

Multicomplex inverse hyperbolic cotangent function via the matrix representation.
"""
function Base.acoth(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = acoth(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    asech(m::Multicomplex)

Multicomplex inverse hyperbolic secant function via the matrix representation.
"""
function Base.asech(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = asech(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end

"""
    acsch(m::Multicomplex)

Multicomplex inverse hyperbolic cosecant function via the matrix representation.
"""
function Base.acsch(m::Multicomplex{T,N,C}) where {T,N,C}
    # Convert to mutable Matrix for LinearAlgebra functions
    mat = Matrix(matrep(m))
    result = acsch(mat)
    Multicomplex{N}(real.(result)[SVector{C}(SOneTo(C))])
end