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