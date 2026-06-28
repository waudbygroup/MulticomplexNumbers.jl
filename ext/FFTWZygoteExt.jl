module FFTWZygoteExt

# The multicomplex FFT operations are implemented with in-place, unsafe
# pointer reinterprets (see `FFTWExt`). Zygote cannot differentiate through that
# implementation, and rather than use the ChainRules `rrule`s it traces the
# in-place code and silently produces identity-like (wrong) gradients.
#
# Defining `Zygote.@adjoint` directly forces Zygote to use the correct adjoint
# transform (the same one expressed by the ChainRules rules in
# `FFTWChainRulesExt`), evaluated with the concrete, working allocating
# transforms:
#
#     fft   ↦  bfft
#     bfft  ↦  fft
#     ifft  ↦  (1/n) * fft        (n = number of transformed points)
#
# Loads only when both `FFTW` and `Zygote` are available.

using MulticomplexNumbers
using MulticomplexNumbers: Multicomplex, order
using FFTW
import FFTW: fft, ifft, bfft
using Zygote: @adjoint

# Number of points transformed along `dims` (integer or collection).
@inline _tlen(A, dims) = prod(size(A, d) for d in (dims isa Integer ? (dims,) : dims))
@inline _u(u::Integer) = u
@inline _u(u::Multicomplex) = order(u)
# Materialise a (possibly lazy, e.g. Fill) cotangent into a dense array, which the
# pointer-based transforms require.
@inline _arr(x::Array) = x
@inline _arr(x::AbstractArray) = collect(x)

# --- fft : adjoint is bfft ---
@adjoint fft(A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}) =
    fft(A, unit), ȳ -> (bfft(_arr(ȳ), _u(unit)), nothing)
@adjoint fft(A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}, dims) =
    fft(A, unit, dims), ȳ -> (bfft(_arr(ȳ), _u(unit), dims), nothing, nothing)

# --- bfft : adjoint is fft ---
@adjoint bfft(A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}) =
    bfft(A, unit), ȳ -> (fft(_arr(ȳ), _u(unit)), nothing)
@adjoint bfft(A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}, dims) =
    bfft(A, unit, dims), ȳ -> (fft(_arr(ȳ), _u(unit), dims), nothing, nothing)

# --- ifft : adjoint is (1/n) * fft ---
@adjoint ifft(A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}) =
    ifft(A, unit), ȳ -> (fft(_arr(ȳ), _u(unit)) ./ length(A), nothing)
@adjoint ifft(A::AbstractArray{<:Multicomplex}, unit::Union{Integer,Multicomplex}, dims) =
    ifft(A, unit, dims), ȳ -> (fft(_arr(ȳ), _u(unit), dims) ./ _tlen(A, dims), nothing, nothing)

end # module
