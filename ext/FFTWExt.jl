module FFTWExt

using MulticomplexNumbers
using FFTW
import FFTW: fft!, ifft!, bfft!, fft, ifft, bfft

"""
    _mc_transform!(transformfn!, A, unit, dims)

Internal helper that performs an FFT-family transform on a multicomplex array.
Converts to complex view, applies `transformfn!`, then converts back.
"""
function _mc_transform!(transformfn!, A::AbstractArray{M}, unit::Integer, dims) where M<:Multicomplex{T,N,C} where {T,N,C}
    w = MulticomplexNumbers.unsafe_ascomplex!(A, unit)
    try
        d = dims .+ (N - 1)
        transformfn!(w, d)
    catch e
        MulticomplexNumbers.unsafe_fromcomplex!(A, unit)
        rethrow(e)
    end
    MulticomplexNumbers.unsafe_fromcomplex!(A, unit)
    return A
end

# --- fft! ---

"""
    fft!(A::AbstractArray{<:Multicomplex}, unit::Integer [, dims])
    fft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex [, dims])

In-place FFT of a multicomplex array, transforming along the specified
imaginary `unit`. The `unit` can be an integer or a multicomplex imaginary
constant (e.g. `im2`), from which `order()` is extracted.

When `dims` is omitted, all array dimensions are transformed.
"""
fft!(A::AbstractArray{<:Multicomplex}, unit::Integer, dims) =
    _mc_transform!(FFTW.fft!, A, unit, dims)

fft!(A::AbstractArray{<:Multicomplex}, unit::Integer) =
    fft!(A, unit, 1:ndims(A))

fft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex, dims) =
    fft!(A, order(unit), dims)

fft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex) =
    fft!(A, order(unit))

# --- ifft! ---

"""
    ifft!(A::AbstractArray{<:Multicomplex}, unit::Integer [, dims])
    ifft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex [, dims])

In-place inverse FFT of a multicomplex array. Same calling conventions as `fft!`.
"""
ifft!(A::AbstractArray{<:Multicomplex}, unit::Integer, dims) =
    _mc_transform!(FFTW.ifft!, A, unit, dims)

ifft!(A::AbstractArray{<:Multicomplex}, unit::Integer) =
    ifft!(A, unit, 1:ndims(A))

ifft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex, dims) =
    ifft!(A, order(unit), dims)

ifft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex) =
    ifft!(A, order(unit))

# --- bfft! ---

"""
    bfft!(A::AbstractArray{<:Multicomplex}, unit::Integer [, dims])
    bfft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex [, dims])

In-place unnormalized inverse FFT (backward FFT). Same calling conventions as `fft!`.
`bfft!(fft!(copy(A), unit), unit) ≈ A .* length(A)`.
"""
bfft!(A::AbstractArray{<:Multicomplex}, unit::Integer, dims) =
    _mc_transform!(FFTW.bfft!, A, unit, dims)

bfft!(A::AbstractArray{<:Multicomplex}, unit::Integer) =
    bfft!(A, unit, 1:ndims(A))

bfft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex, dims) =
    bfft!(A, order(unit), dims)

bfft!(A::AbstractArray{<:Multicomplex}, unit::Multicomplex) =
    bfft!(A, order(unit))

# --- Allocating variants ---

"""
    fft(A::AbstractArray{<:Multicomplex}, unit [, dims])

Allocating FFT — returns a transformed copy, leaving `A` unchanged.
"""
fft(A::AbstractArray{<:Multicomplex}, unit::Integer, dims) =
    fft!(copy(A), unit, dims)

fft(A::AbstractArray{<:Multicomplex}, unit::Integer) =
    fft!(copy(A), unit)

fft(A::AbstractArray{<:Multicomplex}, unit::Multicomplex, dims) =
    fft!(copy(A), order(unit), dims)

fft(A::AbstractArray{<:Multicomplex}, unit::Multicomplex) =
    fft!(copy(A), order(unit))

"""
    ifft(A::AbstractArray{<:Multicomplex}, unit [, dims])

Allocating inverse FFT — returns a transformed copy, leaving `A` unchanged.
"""
ifft(A::AbstractArray{<:Multicomplex}, unit::Integer, dims) =
    ifft!(copy(A), unit, dims)

ifft(A::AbstractArray{<:Multicomplex}, unit::Integer) =
    ifft!(copy(A), unit)

ifft(A::AbstractArray{<:Multicomplex}, unit::Multicomplex, dims) =
    ifft!(copy(A), order(unit), dims)

ifft(A::AbstractArray{<:Multicomplex}, unit::Multicomplex) =
    ifft!(copy(A), order(unit))

"""
    bfft(A::AbstractArray{<:Multicomplex}, unit [, dims])

Allocating unnormalized inverse FFT — returns a transformed copy, leaving `A` unchanged.
"""
bfft(A::AbstractArray{<:Multicomplex}, unit::Integer, dims) =
    bfft!(copy(A), unit, dims)

bfft(A::AbstractArray{<:Multicomplex}, unit::Integer) =
    bfft!(copy(A), unit)

bfft(A::AbstractArray{<:Multicomplex}, unit::Multicomplex, dims) =
    bfft!(copy(A), order(unit), dims)

bfft(A::AbstractArray{<:Multicomplex}, unit::Multicomplex) =
    bfft!(copy(A), order(unit))

end
