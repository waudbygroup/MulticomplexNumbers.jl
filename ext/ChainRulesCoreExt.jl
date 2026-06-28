module ChainRulesCoreExt

# Reverse-mode (and where trivial, forward-mode) differentiation rules for
# multicomplex numbers, enabling Zygote (and any ChainRules-consuming AD) to
# differentiate through multicomplex algebra.
#
# Convention
# ----------
# A `Multicomplex{T,N,C}` is treated as a real vector space of its `C = 2^N`
# components.  The (co)tangent of a multicomplex number is therefore another
# multicomplex number of the same order, holding the partial derivatives with
# respect to each component.  Gradients of real-valued losses with respect to
# multicomplex inputs come back as `Multicomplex` values.
#
# The adjoint of multicomplex multiplication
# ------------------------------------------
# For a real-valued loss the reverse-mode adjoint of "multiply by `y`" is the
# transpose of `y`'s multiplication matrix, `matrep(y)ᵀ`.  One can show
# (and it is verified numerically in the test-suite) that
#
#     matrep(m)ᵀ == matrep(τ(m)),
#
# where `τ` flips the sign of *every* imaginary unit — i.e. it negates each
# component whose index (0-based) has an odd number of set bits.  `τ` is a ring
# automorphism, so for any analytic function `f` (built from `+`, `*`, power
# series, …) the pullback is simply
#
#     m̄ = τ(f′(m)) * z̄ ,   z = f(m).
#
# For `N == 1`, `τ` coincides with ordinary complex conjugation, recovering the
# familiar complex multiplication adjoint `x̄ = z̄ * conj(y)`.

using MulticomplexNumbers
using MulticomplexNumbers: Multicomplex, flat, component, realest, order
using StaticArrays
using ChainRulesCore
using ChainRulesCore: NoTangent, ZeroTangent, AbstractZero, AbstractThunk,
    Tangent, unthunk, ProjectTo, rrule

#####################################################################
# Helpers
#####################################################################

# Sign mask for the "flip every imaginary unit" automorphism τ.
@inline _signmask(::Val{C}) where {C} =
    SVector{C,Int}(ntuple(j -> iseven(count_ones(j - 1)) ? 1 : -1, Val(C)))

"τ(m): negate every component whose index contains an odd number of imaginary units."
@inline _tau(m::Multicomplex{T,N,C}) where {T,N,C} =
    Multicomplex{N}(m.value .* _signmask(Val(C)))

# Canonicalise an incoming cotangent into a `Multicomplex` of matching order.
@inline _asmc(::Multicomplex{T,N,C}, dy::Multicomplex{S,N,C}) where {T,S,N,C} = dy
@inline _asmc(m::Multicomplex{T,N,C}, dy::AbstractThunk) where {T,N,C} = _asmc(m, unthunk(dy))
@inline _asmc(m::Multicomplex{T,N,C}, ::AbstractZero) where {T,N,C} = zero(Multicomplex{T,N,C})
@inline _asmc(m::Multicomplex{T,N,C}, dy::Tangent) where {T,N,C} =
    Multicomplex{N}(SVector{C}(unthunk(getproperty(dy, :value))...))
@inline _asmc(m::Multicomplex{T,N,C}, dy::AbstractVector) where {T,N,C} =
    Multicomplex{N}(SVector{C,T}(dy))
@inline _asmc(::Multicomplex{T,0,1}, dy::Real) where {T} = Multicomplex{0}(SVector(dy))

# One-hot multicomplex with `val` in slot `k`.
@inline _onehot(::Multicomplex{T,N,C}, k::Integer, val) where {T,N,C} =
    Multicomplex{N}(SVector{C}(ntuple(j -> j == k ? val : zero(val), Val(C))))

#####################################################################
# ProjectTo
#####################################################################

ChainRulesCore.ProjectTo(m::Multicomplex{T,N,C}) where {T,N,C} =
    ProjectTo{Multicomplex}(; N = N, C = C)

(p::ProjectTo{Multicomplex})(dx::Multicomplex) = dx
(p::ProjectTo{Multicomplex})(dx::AbstractZero) = dx
(p::ProjectTo{Multicomplex})(dx::AbstractThunk) = p(unthunk(dx))
(p::ProjectTo{Multicomplex})(dx::Tangent) =
    Multicomplex{p.N}(SVector{p.C}(unthunk(getproperty(dx, :value))...))

#####################################################################
# Constructors (entry points from reals / complex / lower-order)
#####################################################################

# Build from a component vector — covers both `Multicomplex{N}(v)` and
# `Multicomplex{T,N,C}(v)` since both are `Type{<:Multicomplex}`.
function rrule(MC::Type{<:Multicomplex}, v::SVector{C,T}) where {C,T}
    y = MC(v)
    construct_pullback(ȳ) = (NoTangent(), _asmc(y, ȳ).value)
    return y, construct_pullback
end

function rrule(::Type{Multicomplex}, x::Real)
    y = Multicomplex(x)
    real0_pullback(ȳ) = (NoTangent(), realest(_asmc(y, ȳ)))
    return y, real0_pullback
end

function rrule(::Type{Multicomplex}, x::Real, y::Real)
    z = Multicomplex(x, y)
    function real1_pullback(z̄)
        c = _asmc(z, z̄)
        return (NoTangent(), c.value[1], c.value[2])
    end
    return z, real1_pullback
end

function rrule(::Type{Multicomplex}, x::Real, y::Real, u::Real, v::Real)
    z = Multicomplex(x, y, u, v)
    function real2_pullback(z̄)
        c = _asmc(z, z̄)
        return (NoTangent(), c.value[1], c.value[2], c.value[3], c.value[4])
    end
    return z, real2_pullback
end

# Build a higher order number from real/imag multicomplex halves.
function rrule(::Type{Multicomplex}, a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    z = Multicomplex(a, b)
    function pair_pullback(z̄)
        c = _asmc(z, z̄)
        ā = Multicomplex{N}(c.value[SOneTo(C)])
        b̄ = Multicomplex{N}(c.value[SOneTo(C) .+ Scalar(C)])
        return (NoTangent(), ā, b̄)
    end
    return z, pair_pullback
end

#####################################################################
# Component access (exit points to reals / lower-order)
#####################################################################

function rrule(::typeof(flat), m::Multicomplex{T,N,C}) where {T,N,C}
    flat_pullback(ȳ) = (NoTangent(), Multicomplex{N}(SVector{C}(unthunk(ȳ)...)))
    return flat(m), flat_pullback
end

function rrule(::typeof(component), m::Multicomplex{T,N,C}, k) where {T,N,C}
    component_pullback(c̄) = (NoTangent(), _onehot(m, k, unthunk(c̄)), NoTangent())
    return component(m, k), component_pullback
end

function rrule(::typeof(realest), m::Multicomplex{T,N,C}) where {T,N,C}
    realest_pullback(c̄) = (NoTangent(), _onehot(m, 1, unthunk(c̄)))
    return realest(m), realest_pullback
end

# real
function rrule(::typeof(real), m::Multicomplex{T,0,1}) where {T}
    real0_pullback(r̄) = (NoTangent(), Multicomplex{0}(SVector(unthunk(r̄))))
    return real(m), real0_pullback
end
function rrule(::typeof(real), m::Multicomplex{T,1,2}) where {T}
    real1_pullback(r̄) = (NoTangent(), Multicomplex{1}(SVector(unthunk(r̄), zero(T))))
    return real(m), real1_pullback
end
function rrule(::typeof(real), m::Multicomplex{T,N,C}) where {T,N,C}
    r = real(m)
    function realN_pullback(r̄)
        cv = _asmc(r, r̄).value
        z = zeros(SVector{C ÷ 2,eltype(cv)})
        return (NoTangent(), Multicomplex{N}(vcat(cv, z)))
    end
    return r, realN_pullback
end

# imag
function rrule(::typeof(imag), m::Multicomplex{T,0,1}) where {T}
    imag0_pullback(ī) = (NoTangent(), ZeroTangent())
    return imag(m), imag0_pullback
end
function rrule(::typeof(imag), m::Multicomplex{T,1,2}) where {T}
    imag1_pullback(ī) = (NoTangent(), Multicomplex{1}(SVector(zero(T), unthunk(ī))))
    return imag(m), imag1_pullback
end
function rrule(::typeof(imag), m::Multicomplex{T,N,C}) where {T,N,C}
    i = imag(m)
    function imagN_pullback(ī)
        cv = _asmc(i, ī).value
        z = zeros(SVector{C ÷ 2,eltype(cv)})
        return (NoTangent(), Multicomplex{N}(vcat(z, cv)))
    end
    return i, imagN_pullback
end

function rrule(::typeof(conj), m::Multicomplex{T,N,C}) where {T,N,C}
    y = conj(m)
    conj_pullback(ȳ) = (NoTangent(), conj(_asmc(y, ȳ)))
    return y, conj_pullback
end

#####################################################################
# Additive structure (linear, trivial adjoints)
#####################################################################

function rrule(::typeof(+), m::Multicomplex{T,N,C}) where {T,N,C}
    y = +m
    uplus_pullback(ȳ) = (NoTangent(), _asmc(y, ȳ))
    return y, uplus_pullback
end
function rrule(::typeof(-), m::Multicomplex{T,N,C}) where {T,N,C}
    y = -m
    uminus_pullback(ȳ) = (NoTangent(), -_asmc(y, ȳ))
    return y, uminus_pullback
end

function rrule(::typeof(+), a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N}
    y = a + b
    function plus_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), c, c)
    end
    return y, plus_pullback
end
function rrule(::typeof(+), a::Multicomplex{T,N}, b::Real) where {T,N}
    y = a + b
    function plusr_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), c, realest(c))
    end
    return y, plusr_pullback
end
function rrule(::typeof(+), a::Real, b::Multicomplex{T,N}) where {T,N}
    y = a + b
    function rplus_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), realest(c), c)
    end
    return y, rplus_pullback
end

function rrule(::typeof(-), a::Multicomplex{T,N}, b::Multicomplex{T,N}) where {T,N}
    y = a - b
    function minus_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), c, -c)
    end
    return y, minus_pullback
end
function rrule(::typeof(-), a::Multicomplex{T,N}, b::Real) where {T,N}
    y = a - b
    function minusr_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), c, -realest(c))
    end
    return y, minusr_pullback
end
function rrule(::typeof(-), a::Real, b::Multicomplex{T,N}) where {T,N}
    y = a - b
    function rminus_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), realest(c), -c)
    end
    return y, rminus_pullback
end

#####################################################################
# Scalar scaling
#####################################################################

function rrule(::typeof(*), s::Real, m::Multicomplex{T,N,C}) where {T,N,C}
    y = s * m
    function smul_pullback(ȳ)
        c = _asmc(y, ȳ)
        s̄ = sum(m.value .* c.value)
        return (NoTangent(), s̄, s * c)
    end
    return y, smul_pullback
end
function rrule(::typeof(*), m::Multicomplex{T,N,C}, s::Real) where {T,N,C}
    y = m * s
    function muls_pullback(ȳ)
        c = _asmc(y, ȳ)
        s̄ = sum(m.value .* c.value)
        return (NoTangent(), s * c, s̄)
    end
    return y, muls_pullback
end
function rrule(::typeof(/), m::Multicomplex{T,N,C}, s::Real) where {T,N,C}
    y = m / s
    function divs_pullback(ȳ)
        c = _asmc(y, ȳ)
        s̄ = -sum(y.value .* c.value) / s
        return (NoTangent(), c / s, s̄)
    end
    return y, divs_pullback
end

#####################################################################
# Multicomplex multiplication, inversion and division
#####################################################################

function rrule(::typeof(*), a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    y = a * b
    function mul_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(b) * c, _tau(a) * c)
    end
    return y, mul_pullback
end

function rrule(::typeof(inv), m::Multicomplex{T,N,C}) where {T,N,C}
    y = inv(m)
    function inv_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(-(y * y)) * c)
    end
    return y, inv_pullback
end

function rrule(::typeof(/), a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
    y = a / b
    ib = inv(b)
    function div_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(ib) * c, _tau(-(y * ib)) * c)
    end
    return y, div_pullback
end

#####################################################################
# Powers
#####################################################################

function rrule(::typeof(^), m::Multicomplex{T,N,C}, p::Integer) where {T,N,C}
    y = m^p
    function powi_pullback(ȳ)
        c = _asmc(y, ȳ)
        m̄ = _tau(p * m^(p - 1)) * c
        return (NoTangent(), m̄, NoTangent())
    end
    return y, powi_pullback
end

function rrule(::typeof(^), m::Multicomplex{T,N,C}, p::Real) where {T,N,C}
    y = m^p
    function powr_pullback(ȳ)
        c = _asmc(y, ȳ)
        m̄ = _tau(p * m^(p - 1)) * c
        p̄ = sum((y * log(m)).value .* c.value)
        return (NoTangent(), m̄, p̄)
    end
    return y, powr_pullback
end

# `m^2`, `m^3`, … are lowered to `Base.literal_pow`.
function rrule(::typeof(Base.literal_pow), ::typeof(^), m::Multicomplex{T,N,C}, ::Val{p}) where {T,N,C,p}
    y = m^p
    function litpow_pullback(ȳ)
        c = _asmc(y, ȳ)
        m̄ = _tau(p * m^(p - 1)) * c
        return (NoTangent(), NoTangent(), m̄, NoTangent())
    end
    return y, litpow_pullback
end

#####################################################################
# Absolute values
#####################################################################

function rrule(::typeof(abs2), m::Multicomplex{T,N,C}) where {T,N,C}
    y = abs2(m)
    function abs2_pullback(ȳ)
        s̄ = unthunk(ȳ)
        return (NoTangent(), Multicomplex{N}((2 * s̄) .* m.value))
    end
    return y, abs2_pullback
end

function rrule(::typeof(abs), m::Multicomplex{T,N,C}) where {T,N,C}
    y = abs(m)
    function abs_pullback(ȳ)
        s̄ = unthunk(ȳ)
        return (NoTangent(), Multicomplex{N}((s̄ / y) .* m.value))
    end
    return y, abs_pullback
end

#####################################################################
# Elementary transcendental functions: m̄ = τ(f′(m)) * z̄
#####################################################################

for (f, dfexpr) in (
    (:exp,  :(y)),               # exp′ = exp
    (:sin,  :(cos(m))),
    (:cos,  :(-sin(m))),
    (:sinh, :(cosh(m))),
    (:cosh, :(sinh(m))),
)
    @eval function rrule(::typeof($f), m::Multicomplex{T,N,C}) where {T,N,C}
        y = $f(m)
        function $(Symbol(f, :_pullback))(ȳ)
            c = _asmc(y, ȳ)
            return (NoTangent(), _tau($dfexpr) * c)
        end
        return y, $(Symbol(f, :_pullback))
    end
end

function rrule(::typeof(log), m::Multicomplex{T,N,C}) where {T,N,C}
    y = log(m)
    function log_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(inv(m)) * c)
    end
    return y, log_pullback
end

function rrule(::typeof(sqrt), m::Multicomplex{T,N,C}) where {T,N,C}
    y = sqrt(m)
    function sqrt_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(inv(2 * y)) * c)
    end
    return y, sqrt_pullback
end

function rrule(::typeof(tan), m::Multicomplex{T,N,C}) where {T,N,C}
    y = tan(m)
    function tan_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(one(m) + y * y) * c)   # sec² = 1 + tan²
    end
    return y, tan_pullback
end

function rrule(::typeof(tanh), m::Multicomplex{T,N,C}) where {T,N,C}
    y = tanh(m)
    function tanh_pullback(ȳ)
        c = _asmc(y, ȳ)
        return (NoTangent(), _tau(one(m) - y * y) * c)   # sech² = 1 - tanh²
    end
    return y, tanh_pullback
end

end # module
