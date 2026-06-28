module ZygoteChainRulesExt

# Zygote ships its own native adjoints for a handful of functions
# (`real`, `imag`, `inv`, `^`, `Base.literal_pow`) that treat any `Number` —
# including a `Multicomplex` — as if it were a complex scalar.  Those native
# adjoints take precedence over the ChainRules `rrule`s and are therefore wrong
# for multicomplex numbers (wrong sign on the imaginary components, or a
# truncated-size gradient).
#
# This extension overrides them for `Multicomplex` arguments by delegating to
# the ChainRules `rrule` defined in `ChainRulesCoreExt`, keeping a single source
# of truth for the differentiation rules.  It loads only when both
# `ChainRulesCore` and `Zygote` are available.

using MulticomplexNumbers
using MulticomplexNumbers: Multicomplex
using Zygote: @adjoint
import ChainRulesCore
using ChainRulesCore: AbstractZero

# Zygote expects `nothing` (not a ChainRules zero) for a non-differentiable slot.
@inline _z(x) = x
@inline _z(::AbstractZero) = nothing

@adjoint function Base.inv(m::Multicomplex)
    y, back = ChainRulesCore.rrule(inv, m)
    return y, ȳ -> (_z(back(ȳ)[2]),)
end

@adjoint function Base.real(m::Multicomplex)
    y, back = ChainRulesCore.rrule(real, m)
    return y, ȳ -> (_z(back(ȳ)[2]),)
end

@adjoint function Base.imag(m::Multicomplex)
    y, back = ChainRulesCore.rrule(imag, m)
    return y, ȳ -> (_z(back(ȳ)[2]),)
end

@adjoint function Base.:^(m::Multicomplex, p::Real)
    y, back = ChainRulesCore.rrule(^, m, p)
    return y, ȳ -> begin
        b = back(ȳ)
        (_z(b[2]), _z(b[3]))
    end
end

@adjoint function Base.literal_pow(::typeof(^), m::Multicomplex, v::Val)
    y, back = ChainRulesCore.rrule(Base.literal_pow, ^, m, v)
    return y, ȳ -> begin
        b = back(ȳ)
        (nothing, _z(b[3]), nothing)
    end
end

end # module
