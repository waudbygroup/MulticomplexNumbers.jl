module FFTWExt

"""
notes on testing:
]activate --temp
]dev .
]add FFTW
using MulticomplexNumbers, FFTW
"""

using MulticomplexNumbers
using FFTW
import FFTW.fft!
@info "FFTWExt being loaded"



"""
outline of approach:

fft - need to copy array, do in-place fft!, re-interpret again without copying
fft! - don't copy array, do in-place fft
"""

function fft!(A::AbstractArray{M}, unit::Integer) where M<:Multicomplex{T,N,C} where {T,N,C}
    fft!(A, unit, 1:ndims(A))
end

function fft!(A::AbstractArray{M}, unit::Integer, dims) where M<:Multicomplex{T,N,C} where {T,N,C}
    w = MulticomplexNumbers.unsafe_ascomplex!(A, unit)
    try
        if N == 1
            fft!(w, dims)
        elseif N == 2
            d = dims .+ 1
            fft!(w, d)
        elseif N == 3
            d = dims .+ 2
            fft!(w, d)
        end
    catch e
        MulticomplexNumbers.unsafe_fromcomplex!(A, unit)
        throw(e)
    end
    MulticomplexNumbers.unsafe_fromcomplex!(A, unit)
end

end