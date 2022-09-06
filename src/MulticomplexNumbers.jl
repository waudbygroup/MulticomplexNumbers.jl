# module MulticomplexNumbers

using StaticArrays

# export greet

# """Greetings!"""
# greet() = println("Hello world!")


################
# Multicomplex #
################

"""
    can_multicomplex(::Type)
Determines whether a type is allowed as the scalar type in a
multicomplex number. By default, only `<:Real` types are allowed.
"""
can_multicomplex(::Type{<:Real}) = true
can_multicomplex(::Type) = false


"""
    Multicomplex{T,N,C}
Defines a multicomplex number ℂ_N, over base field T, having C=2^N components.
"""
struct Multicomplex{T,N,C} <: Number
    value::SVector{C,T}
    function Multicomplex{N}(value::SVector{C, T}) where {T, N, C}
        can_multicomplex(T) || throw_cannot_multicomplex(T)
        2^N == C || throw(ArgumentError("Cannot create a multicomplex C_$N with $C components (expect C=2^N=$(2^N)) components)"))
        new{T, N, C}(value)
    end
end

"Build a higher order multicomplex number from two multicomplex 'real' and 'imaginary' components"
function Multicomplex(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N+1}(vcat(a.value,b.value))
end

# Multicomplex{0} constructors
Multicomplex(x::Real) = Multicomplex{0}(SVector(x))

# Multicomplex{1} constructors
Multicomplex(x::Real, y::Real) = Multicomplex{1}(SVector(x,y))
Multicomplex(z::Complex) = Multicomplex{1}(SVector(z.re, z.im))

# Multicomplex{2} constructors
"Multicomplex(x,y,u,v) = x + i¹y + i²u + i¹i²v"
Multicomplex(x::Real, y::Real, u::Real, v::Real) = Multicomplex{2}(SVector(promote(x,y,u,v)...))
Multicomplex(a::Complex, b::Complex) = Multicomplex{2}(SVector(a.re, a.im, b.re, b.im))
Multicomplex(a::Complex, b::Real) = Multicomplex(a, complex(b))
Multicomplex(a::Real, b::Complex) = Multicomplex(complex(a), b)


"""
    im1, im2, im3, im4, im5, im6
The imaginary units, up to dimension 6.
"""
const im1 = Multicomplex{1}(SVector{2,Int8}(0, 1))
const im2 = Multicomplex{2}(SVector{4,Int8}(0, 0, 1, 0))
const im3 = Multicomplex{3}(SVector{8,Int8}(0, 0, 0, 0, 1, 0, 0, 0))
const im4 = Multicomplex{4}(SVector{16,Int8}(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
const im5 = Multicomplex{5}(SVector{32,Int8}(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
const im6 = Multicomplex{6}(SVector{64,Int8}(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))


############################
# Conversion and promotion #
############################

# cast real into defined multicomplex type:
Multicomplex{T}(x::Real) where {T<:Real} = Multicomplex(convert(T,x))
# cast complex into defined multicomplex type:
Multicomplex{T}(z::Complex) where {T<:Real} = Multicomplex(reim(convert(Complex{T}, z))...)
# change multicomplex type:
Multicomplex{T,N,C}(m::Multicomplex{S,N,C}) where {T<:Real,N,C,S} = Multicomplex{N}(convert(SVector{C,T}, m.value))
# convert multicomplex size and/or type:
Multicomplex{T,1,2}(m::Multicomplex{S,0,1}) where {T,S} = Multicomplex{1}(SVector{2,T}(m.value[1], zero(T)))
Multicomplex{T,2,4}(m::Multicomplex{S,0,1}) where {T,S} = Multicomplex{2}(SVector{4,T}(m.value[1], zero(T), zero(T), zero(T)))
Multicomplex{T,2,4}(m::Multicomplex{S,1,2}) where {T,S} = Multicomplex{2}(SVector{4,T}(m.value[1], m.value[2], zero(T), zero(T)))
Multicomplex{T,3,8}(m::Multicomplex{S,0,1}) where {T,S} = Multicomplex{3}(SVector{8,T}(m.value[1], zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T)))
Multicomplex{T,3,8}(m::Multicomplex{S,1,2}) where {T,S} = Multicomplex{3}(SVector{8,T}(m.value[1], m.value[2], zero(T), zero(T), zero(T), zero(T), zero(T), zero(T)))
Multicomplex{T,3,8}(m::Multicomplex{S,2,4}) where {T,S} = Multicomplex{3}(SVector{8,T}(m.value[1], m.value[2], m.value[3], m.value[4], zero(T), zero(T), zero(T), zero(T)))
function Multicomplex{T,N,C}(m::Multicomplex{S,N0,C0}) where {T<:Real,N,C,S,N0,C0}
    N0 ≤ N || throw(InexactError(nameof(Multicomplex),Multicomplex{T,N,C},m))
    v = zeros(MVector{C,T})
    for i=1:C0
        v[i] = m.value[i]
    end
    Multicomplex{N}(SVector(v))
end

(::Type{T})(m::Multicomplex) where {T<:Real} = isreal(m) ? T(realest(m))::T : throw(InexactError(nameof(T), T, m))

Multicomplex(m::Multicomplex) = m

Base.promote_rule(::Type{Multicomplex{T,N,C}}, ::Type{S}) where {T<:Real,S<:Real,N,C} =
    Multicomplex{promote_type(T,S),N,C}
Base.promote_rule(::Type{Multicomplex{T,N,C}}, ::Type{Complex{S}}) where {T<:Real,S<:Real,N,C} =
    Multicomplex{promote_type(T,S),N,C}
Base.promote_rule(
        ::Type{<:Multicomplex{T,N1,C1}},
        ::Type{<:Multicomplex{S,N2,C2}}) where
        {T<:Real,S<:Real,N1,C1,N2,C2} =
    Multicomplex{promote_type(T,S),max(N1,N2),max(C1,C2)}

Base.widen(::Type{Multicomplex{T,N,C}}) where {T,N,C} = Multicomplex{widen(T),N,C}
Base.float(::Type{Multicomplex{T,N,C}}) where {T<:AbstractFloat,N,C} = Multicomplex{T,N,C}
Base.float(::Type{Multicomplex{T,N,C}}) where {T,N,C} = Multicomplex{float(T),N,C}



##############
# Components #
##############

"""Utility function to extract a real-valued component from a multicomplex number"""
component(m::Multicomplex, k) = m.value[k]

"""Extract the 'most real' component"""
realest(m::Multicomplex) = m.value[1]

Base.real(m::Multicomplex{T,0}) where {T} = m.value[1]
Base.real(m::Multicomplex{T,1}) where {T} = m.value[1]
Base.real(m::Multicomplex{T,N,C}) where {T,N,C} = Multicomplex{N-1}(m.value[SOneTo(C÷2)])

Base.imag(::Multicomplex{T,0}) where {T} = zero(T)
Base.imag(m::Multicomplex{T,1}) where {T} = m.value[2]
Base.imag(m::Multicomplex{T,N,C}) where {T,N,C} = Multicomplex{N-1}(m.value[SOneTo(C÷2) .+ Scalar(C÷2)])

"""testing for realness"""
Base.isreal(m::Multicomplex{T,N,C}) where {T,N,C} = iszero(m.value[SOneTo(C-1) .+ Scalar(1)])
Base.isinteger(m::Multicomplex) = isreal(m) & isinteger(m.value[1])
Base.isfinite(m::Multicomplex) = isfinite(m.value)
Base.isnan(m::Multicomplex) = isnan(m.value)
Base.isinf(m::Multicomplex) = isinf(m.value)
Base.iszero(m::Multicomplex) = iszero(m.value)
Base.isone(m::Multicomplex) = isreal(m) & isone(m.value[1])


##########################
# Matrix representations #
##########################

"""matrix representation of a multicomplex number"""
matrep(m::Multicomplex{T,0}) where {T} = SMatrix{1,1,T}(m.value[1])
matrep(m::Multicomplex{T,1}) where {T} = SMatrix{2,2,T}(m.value[1], m.value[2], -m.value[2], m.value[1])
matrep(m::Multicomplex{T,2}) where {T} = SMatrix{4,4,T}(m.value[1], m.value[2], m.value[3], m.value[4], -m.value[2], m.value[1], -m.value[4], m.value[3], -m.value[3], -m.value[4], m.value[1], m.value[2], m.value[4], -m.value[3], -m.value[2], m.value[1])
    # m.value[1], -m.value[2], -m.value[3],  m.value[4],  m.value[2],  m.value[1], -m.value[4], -m.value[3],  m.value[3], -m.value[4],  m.value[1], -m.value[2],  m.value[4],  m.value[3],  m.value[2],  m.value[1])
matrep(m::Multicomplex{T,3}) where {T} = SMatrix{8,8,T}(m.value[1], m.value[2], m.value[3], m.value[4], m.value[5], m.value[6], m.value[7], m.value[8], -m.value[2], m.value[1], -m.value[4], m.value[3], -m.value[6], m.value[5], -m.value[8], m.value[7], -m.value[3], -m.value[4], m.value[1], m.value[2], -m.value[7], -m.value[8], m.value[5], m.value[6], m.value[4], -m.value[3], -m.value[2], m.value[1], m.value[8], -m.value[7], -m.value[6], m.value[5], -m.value[5], -m.value[6], -m.value[7], -m.value[8], m.value[1], m.value[2], m.value[3], m.value[4], m.value[6], -m.value[5], m.value[8], -m.value[7], -m.value[2], m.value[1], -m.value[4], m.value[3], m.value[7], m.value[8], -m.value[5], -m.value[6], -m.value[3], -m.value[4], m.value[1], m.value[2], -m.value[8], m.value[7], m.value[6], -m.value[5], m.value[4], -m.value[3], -m.value[2], m.value[1])
# SMatrix{8,8,T}(m.value[1] , -m.value[2] , -m.value[3] , m.value[4] , -m.value[5] , m.value[6] , m.value[7] , -m.value[8] , m.value[2] , m.value[1] , -m.value[4] , -m.value[3] , -m.value[6] , -m.value[5] , m.value[8] , m.value[7] , m.value[3] , -m.value[4] , m.value[1] , -m.value[2] , -m.value[7] , m.value[8] , -m.value[5] , m.value[6] , m.value[4] , m.value[3] , m.value[2] , m.value[1] , -m.value[8] , -m.value[7] , -m.value[6] , -m.value[5] , m.value[5] , -m.value[6] , -m.value[7] , m.value[8] , m.value[1] , -m.value[2] , -m.value[3] , m.value[4] , m.value[6] , m.value[5] , -m.value[8] , -m.value[7] , m.value[2] , m.value[1] , -m.value[4] , -m.value[3] , m.value[7] , -m.value[8] , m.value[5] , -m.value[6] , m.value[3] , -m.value[4] , m.value[1] , -m.value[2] , m.value[8] , m.value[7] , m.value[6] , m.value[5] , m.value[4] , m.value[3] , m.value[2] , m.value[1])
function matrep(m::Multicomplex{T,N,C}) where {T,N,C}
    A = MMatrix{C,C,T}(undef)
    r = matrep(real(m))
    i = matrep(imag(m))
    idx1 = SOneTo(C÷2)
    idx2 = SOneTo(C÷2) .+ Scalar(C÷2)
    A[idx1, idx1] .= r
    A[idx1, idx2] .= -i
    A[idx2, idx1] .= i
    A[idx2, idx2] .= r
    SMatrix(A)
end


"""
    ascomplex(A::AbstractArray{M})

Returns a view of the multicomplex input array A as an array of complex numbers, mapping i⁽¹⁾ -> im

If A has size (m, n, ...) with multicomplex numbers of dimension N, then the output
will have size (2^(N-1), m, n, ...).
"""
function ascomplex(A::AbstractArray{M}) where M<:Multicomplex{T} where {T}
    reinterpret(reshape,Complex{T}, A)
end



##############
# Arithmetic #
##############

# invert if y<0
Base.flipsign(x::Multicomplex, y::Real) = ifelse(signbit(y), -x, x)

# byte order swaps: components are swapped individually
Base.bswap(m::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(bswap.(m.value))

Base.:(==)(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = (a.value == b.value)
# Base.:(==)(a::Multicomplex{T,N1,C1}, b::Multicomplex{S,N2,C2}) where {T,S,N1,N2,C1,C2} = ==(promote(a,b)...)

Base.isequal(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = isequal(a.value, b.value)
# Base.isequal(a::Multicomplex{T,N1,C1}, b::Multicomplex{S,N2,C2}) where {T,S,N1,N2,C1,C2} = isequal(promote(a,b)...)

Base.hash(m::Multicomplex, h::UInt) = hash(m.value, h)

# addition and subtraction
Base.:(+)(a::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(a.value)
Base.:(-)(a::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(-a.value)
Base.:(+)(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(a.value + b.value)
Base.:(-)(a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(a.value - b.value)

# multiplication by a scalar
Base.:(*)(a::Real, m::Multicomplex{T,N,C}) where {T,N,C} = Multicomplex{N}(a*m.value)
Base.:(*)(m::Multicomplex{T,N,C}, a::Real) where {T,N,C} = Multicomplex{N}(a*m.value)

# multicomplex multiplication via matrix representation
function Base.:(*)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    Multicomplex{N}((matrep(a) * matrep(b))[SVector{C}(SOneTo(C))])
    # p = (matrep(a) * matrep(b))
    # v = p[SVector{C}(SOneTo(C))]
    # Multicomplex{N}(v)
end

# conjugation
conj(m::Multicomplex) = Multicomplex(real(m),-imag(m))

##########
# Traits #
##########
Base.ArithmeticStyle(::Type{<:Multicomplex{T}}) where {T} = Base.ArithmeticStyle(T)





##############
# Exceptions #
##############

@noinline function throw_cannot_multicomplex(T::Type)
    throw(ArgumentError("Cannot create a multicomplex number over scalar type $T." *
        " If the type behaves as a scalar, define MulticomplexNumbers.can_multicomplex(::Type{$T}) = true."))
end

# end
