using Random

"""
    can_multicomplex(::Type)
Determines whether a type is allowed as the scalar type in a
multicomplex number. By default, only `<:Real` types are allowed.
"""
can_multicomplex(::Type{<:Real}) = true
can_multicomplex(::Type) = false


"""
    Multicomplex{T,N,C}(value)
Defines a multicomplex number ℂ_N, over base field T, having C=2^N components specified in static vector `value`.

Input:
    value::SVector{C, T})
"""
struct Multicomplex{T,N,C} <: Number
    value::SVector{C,T}
    function Multicomplex{T,N,C}(value::SVector{C, T}) where {T, N, C}
        can_multicomplex(T) || throw_cannot_multicomplex(T)
        2^N == C || throw(ArgumentError("Cannot create a multicomplex C_$N with $C components (expect C=2^N=$(2^N)) components)"))
        new{T, N, C}(value)
    end
end

"""
    Multicomplex{N}(value)
Defines a multicomplex number ℂ_N, with base field and components defined in value.

Input:
    value::SVector{C, T})
"""
Multicomplex{N}(value::SVector{C, T}) where {T, N, C} = Multicomplex{T, N, C}(value)

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
Multicomplex(x::Real, y::Real, u::Real, v::Real) = Multicomplex{2}(SVector(x,y,u,v))
Multicomplex(a::Complex, b::Complex) = Multicomplex{2}(SVector(a.re, a.im, b.re, b.im))
Multicomplex(a::Complex, b::Real) = Multicomplex(a, complex(b))
Multicomplex(a::Real, b::Complex) = Multicomplex(complex(a), b)

# Multicomplex{3} constructor from reals
Multicomplex(rrr::Real, irr::Real, rir::Real, iir::Real, rri::Real, iri::Real, rii::Real, iii::Real) = Multicomplex{3}(SVector(rrr, irr, rir, iir, rri, iri, rii, iii))

# Multicomplex{4} constructor from reals
Multicomplex(rrrr::Real, irrr::Real, rirr::Real, iirr::Real, rrir::Real, irir::Real, riir::Real, iiir::Real,
             rrri::Real, irri::Real, riri::Real, iiri::Real, rrii::Real, irii::Real, riii::Real, iiii::Real) = Multicomplex{4}(SVector(rrrr, irrr, rirr, iirr, rrir, irir, riir, iiir, rrri, irri, riri, iiri, rrii, irii, riii, iiii))


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
# Multicomplex{T,N,C}(x::Real) where {T<:Real,N,C} = Multicomplex(convert(T,x))
Multicomplex{T,0,1}(x::Real) where {T<:Real} = Multicomplex(convert(T,x))
Multicomplex{T,N,C}(x::Real) where {T<:Real,N,C} = Multicomplex{T,N,C}(convert(Multicomplex{T,0,1},x))

# cast complex into defined multicomplex type:
Multicomplex{T,1,2}(z::Complex) where {T<:Real} = Multicomplex(reim(convert(Complex{T}, z))...)
Multicomplex{T,N,C}(z::Complex) where {T<:Real,N,C} = Multicomplex{T,N,C}(convert(Multicomplex{T,1,2},z))

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

Base.rtoldefault(::Type{Multicomplex{T,N,C}}) where {T,N,C} = Base.rtoldefault(T)
Base.rtoldefault(::Multicomplex{T,N,C}) where {T,N,C} = Base.rtoldefault(T)
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:Multicomplex,S<:Number}
    rtol = max(Base.rtoldefault(T),Base.rtoldefault(real(S)))
    return atol > 0 ? zero(rtol) : rtol
end
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:Number,S<:Multicomplex}
    rtol = max(Base.rtoldefault(real(T)),Base.rtoldefault(S))
    return atol > 0 ? zero(rtol) : rtol
end
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:Multicomplex,S<:Multicomplex}
    rtol = max(Base.rtoldefault(T),Base.rtoldefault(S))
    return atol > 0 ? zero(rtol) : rtol
end

##############
# Components #
##############

"""Order of the multicomplex number is the number of imaginary units"""
order(::Real) = 0
order(::Complex) = 1
order(::Multicomplex{T,N}) where {T,N} = N

"""Get components as a vector"""
flat(x::Real) = SVector(x)
flat(z::Complex) = SVector(real(z), imag(z))
flat(m::Multicomplex) = m.value

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
Base.isinteger(m::Multicomplex) = isreal(m) && isinteger(m.value[1])
Base.isfinite(m::Multicomplex) = all(isfinite.(m.value))
Base.isnan(m::Multicomplex) = any(isnan.(m.value))
Base.isinf(m::Multicomplex) = any(isinf.(m.value))
Base.iszero(m::Multicomplex) = all(iszero.(m.value))
Base.isone(m::Multicomplex) = isreal(m) && isone(m.value[1])


"""Norm is sum of squares of components"""
LinearAlgebra.norm(m::Multicomplex) = norm(m.value)


##########
# Traits #
##########
Base.ArithmeticStyle(::Type{<:Multicomplex{T}}) where {T} = Base.ArithmeticStyle(T)


################
# Constructors #
################

"""
    zero(::Type{Multicomplex{T,N,C}})

Return the additive identity (zero) for the given multicomplex type.
"""
Base.zero(::Type{Multicomplex{T,N,C}}) where {T,N,C} = Multicomplex{N}(zeros(SVector{C,T}))
Base.zero(m::Multicomplex{T,N,C}) where {T,N,C} = zero(Multicomplex{T,N,C})

"""
    one(::Type{Multicomplex{T,N,C}})

Return the multiplicative identity (one) for the given multicomplex type.
"""
function Base.one(::Type{Multicomplex{T,N,C}}) where {T,N,C}
    v = zeros(MVector{C,T})
    v[1] = one(T)
    Multicomplex{N}(SVector(v))
end
Base.one(m::Multicomplex{T,N,C}) where {T,N,C} = one(Multicomplex{T,N,C})


######################
# Random number generation
######################

"""
    rand([rng=GLOBAL_RNG], ::Type{Multicomplex{T,N,C}})

Generate a random multicomplex number with components drawn from the default distribution for type T.

# Examples
```julia
julia> rand(Multicomplex{Float64,1,2})
0.234 + 0.891*im1

julia> rand(Multicomplex{Float64,2,4})
(0.123 + 0.456*im1) + (0.789 + 0.234*im1)*im2
```
"""
function Base.rand(rng::AbstractRNG, ::Random.SamplerType{Multicomplex{T,N,C}}) where {T,N,C}
    Multicomplex{N}(rand(rng, SVector{C,T}))
end

"""
    randn([rng=GLOBAL_RNG], ::Type{Multicomplex{T,N,C}}) where {T<:AbstractFloat}

Generate a random multicomplex number with components drawn from a standard normal distribution.

Only available for floating-point types.

# Examples
```julia
julia> randn(Multicomplex{Float64,1,2})
-0.543 + 1.234*im1

julia> randn(Multicomplex{Float64,2,4})
(0.891 - 0.234*im1) + (-1.567 + 0.432*im1)*im2
```
"""
function Base.randn(rng::AbstractRNG, ::Type{Multicomplex{T,N,C}}) where {T<:AbstractFloat,N,C}
    Multicomplex{N}(SVector{C,T}(randn(rng, T) for _ in 1:C))
end

# Convenience methods without explicit RNG
Base.randn(::Type{Multicomplex{T,N,C}}) where {T<:AbstractFloat,N,C} = randn(Random.default_rng(), Multicomplex{T,N,C})


##############
# Exceptions #
##############

@noinline function throw_cannot_multicomplex(T::Type)
    throw(ArgumentError("Cannot create a multicomplex number over scalar type $T." *
        " If the type behaves as a scalar, define MulticomplexNumbers.can_multicomplex(::Type{$T}) = true."))
end