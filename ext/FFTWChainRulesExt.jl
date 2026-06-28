module FFTWChainRulesExt

# Reverse-mode differentiation rules for the multicomplex FFT operations defined
# in `FFTWExt`.  This extension loads only when both `FFTW` and `ChainRulesCore`
# are available.
#
# The transforms are linear, so each pullback is the adjoint transform along the
# same imaginary unit and dimensions (verified numerically in the test-suite):
#
#     fft   ↦  bfft
#     bfft  ↦  fft
#     ifft  ↦  (1/n) * fft        (n = number of transformed points)
#
# Only the allocating variants (`fft`, `ifft`, `bfft`) are differentiable; the
# in-place `fft!`/`ifft!`/`bfft!` mutate their argument and are unsupported by
# reverse-mode AD.

using MulticomplexNumbers
using MulticomplexNumbers: Multicomplex, order
using FFTW
import FFTW: fft, ifft, bfft
using ChainRulesCore
using ChainRulesCore: NoTangent, unthunk, rrule

# Number of points transformed along `dims` (which may be an integer or a collection).
@inline _tlen(A, dims) = prod(size(A, d) for d in (dims isa Integer ? (dims,) : dims))

_unit(u::Integer) = u
_unit(u::Multicomplex) = order(u)

# --- fft : pullback is bfft ---

function rrule(::typeof(fft), A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}, dims)
    y = fft(A, unit, dims)
    u = _unit(unit)
    fft_pullback(ȳ) = (NoTangent(), bfft(unthunk(ȳ), u, dims), NoTangent(), NoTangent())
    return y, fft_pullback
end

function rrule(::typeof(fft), A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex})
    y = fft(A, unit)
    u = _unit(unit)
    fft_pullback(ȳ) = (NoTangent(), bfft(unthunk(ȳ), u), NoTangent())
    return y, fft_pullback
end

# --- bfft : pullback is fft ---

function rrule(::typeof(bfft), A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}, dims)
    y = bfft(A, unit, dims)
    u = _unit(unit)
    bfft_pullback(ȳ) = (NoTangent(), fft(unthunk(ȳ), u, dims), NoTangent(), NoTangent())
    return y, bfft_pullback
end

function rrule(::typeof(bfft), A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex})
    y = bfft(A, unit)
    u = _unit(unit)
    bfft_pullback(ȳ) = (NoTangent(), fft(unthunk(ȳ), u), NoTangent())
    return y, bfft_pullback
end

# --- ifft : pullback is (1/n) * fft ---

function rrule(::typeof(ifft), A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}, dims)
    y = ifft(A, unit, dims)
    u = _unit(unit)
    n = _tlen(A, dims)
    ifft_pullback(ȳ) = (NoTangent(), fft(unthunk(ȳ), u, dims) ./ n, NoTangent(), NoTangent())
    return y, ifft_pullback
end

function rrule(::typeof(ifft), A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex})
    y = ifft(A, unit)
    u = _unit(unit)
    n = length(A)
    ifft_pullback(ȳ) = (NoTangent(), fft(unthunk(ȳ), u) ./ n, NoTangent())
    return y, ifft_pullback
end

end # module
