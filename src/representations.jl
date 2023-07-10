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
    ascomplex(A::AbstractArray{M}, [unit])

Returns a view of the multicomplex input array A as an array of complex numbers, mapping i(unit) -> im.
If a unit is not supplied, defaults to the highest order.

If A has size (m, n, ...) with multicomplex numbers of dimension N, then the output
will have size (2^(N-1), m, n, ...).
"""
function ascomplex(A::AbstractArray{M}, unit::Integer) where M<:Multicomplex{T,N,C} where {T,N,C}
    unit > 0 || throw(ArgumentError("unit must be positive"))
    unit <= N || throw(ArgumentError("unit must be less than multicomplex order"))

    if N == 1
        return reinterpret(reshape, Complex{T}, A)
    elseif N == 2
        if unit == 1
            reinterpret(reshape, Complex{T}, A)
        else
            tmp = reinterpret(reshape, T, A)
            tmp = reshape(tmp, (2,2,size(A)...))
            tmp = permutedims(tmp, (2,1,(2 .+ (1:ndims(A)))...))
            reinterpret(reshape, Complex{T}, tmp)
        end
    elseif N == 3
        if unit == 1
            reinterpret(reshape, Complex{T}, A)
        elseif unit == 2
            tmp = reinterpret(reshape, T, A)
            tmp = reshape(tmp, (2,2,2,size(A)...))
            tmp = permutedims(tmp, (2,1,3,(3 .+ (1:ndims(A)))...))
            tmp = reinterpret(reshape, Complex{T}, tmp)
            reshape(tmp, (4,size(A)...))
        else
            tmp = reinterpret(reshape, T, A)
            tmp = reshape(tmp, (2,2,2,size(A)...))
            tmp = permutedims(tmp, (3,1,2,(3 .+ (1:ndims(A)))...))
            reinterpret(reshape, Complex{T}, tmp)
            tmp = reinterpret(reshape, Complex{T}, tmp)
            reshape(tmp, (4,size(A)...))
        end
    else
        if unit == 1
            reinterpret(reshape, Complex{T}, A)
        else
            tmp = reinterpret(reshape, T, A)
            d = repeat([2],N)
            tmp = reshape(tmp, (d...,size(A)...))
            # unit = 2  ==>  2, 1, 3, 4...
            # unit = 3  ==>  3, 1, 2, 4...
            u = collect(1:N)
            popat!(u, unit)
            pushfirst!(u, unit)
            tmp = permutedims(tmp, (u...,(N .+ (1:ndims(A)))...))
            tmp = reinterpret(reshape, Complex{T}, tmp)
            reshape(tmp, (2^(N-1),size(A)...))
        end
    end
end


function ascomplex(A::AbstractArray{M}) where M<:Multicomplex{T,N,C} where {T,N,C}
    ascomplex(A, N)
end


function unsafe_ascomplex!(A::AbstractArray{M}, unit::Integer) where M<:Multicomplex{T,N,C} where {T,N,C}
    unit > 0 || throw(ArgumentError("unit must be positive"))
    unit <= N || throw(ArgumentError("unit must be less than multicomplex order"))

    if N == 1
        _unsafe_ascomplex_11!(A)
    elseif N == 2
        if unit == 1
            _unsafe_ascomplex_21!(A)
        else
            _unsafe_ascomplex_22!(A)
        end
    elseif N == 3
        if unit == 1
            _unsafe_ascomplex_31!(A)
        elseif unit == 2
            _unsafe_ascomplex_32!(A)
        else
            _unsafe_ascomplex_33!(A)
        end
    else
        throw(ArgumentError("unsupported multicomplex order"))
    end
end

"""
    _unsafe_ascomplex_11(A)

N = 1: (r, i)
im1: (r, i) -> (complex)
"""
function _unsafe_ascomplex_11!(A::AbstractArray{M}) where M<:Multicomplex{T,1} where {T}
    p = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p, size(A))
end


"""
    _unsafe_ascomplex_21(A)

N = 2: (rr, ir, ri, ii)
im1: (rr, ir, ri, ii) -> (complex, 2) [zr, zi]    
"""
function _unsafe_ascomplex_21!(A::AbstractArray{M}) where M<:Multicomplex{T,2} where {T}
    p = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p, (2, size(A)...))
end

"""
    _unsafe_ascomplex_22(A)

N = 2: (rr, ir, ri, ii)
im2: (rr, ir, ri, ii) -> permute 2/3 -> (rr, ri, ir, ii) -> (complex, 2) [zr, zi]
"""
function _unsafe_ascomplex_22!(A::AbstractArray{M}) where M<:Multicomplex{T,2} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (4, length(A)))  # 4 x N
    @inbounds for k = 1:length(A)
        w[2,k], w[3,k] = w[3,k], w[2,k]
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, size(A)...))
end



"""
    _unsafe_ascomplex_31(A)

N = 3: (rrr, irr, rir, iir, rri, iri, rii, iii)
im1: (rrr, irr, rir, iir, rri, iri, rii, iii) -> (complex, 2, 2) [zrr, zir, zri, zii]
"""
function _unsafe_ascomplex_31!(A::AbstractArray{M}) where M<:Multicomplex{T,3} where {T}
    p = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p, (2, 2, size(A)...))
end

"""
    _unsafe_ascomplex_32(A)

N = 3: (rrr, irr, rir, iir, rri, iri, rii, iii)
im2: (rrr, irr, rir, iir, rri, iri, rii, iii) -> permute 2/3, 6/7
    -> (rrr, rir, irr, iir, rri, rii, iri, iii) -> (complex, 2, 2) [zrr, zir, zri, zii]
"""
function _unsafe_ascomplex_32!(A::AbstractArray{M}) where M<:Multicomplex{T,3} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (4, 2*length(A)))  # 4 x N
    @inbounds for k = 1:(2*length(A))
        w[2,k], w[3,k] = w[3,k], w[2,k]
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, 2, size(A)...))
end

"""
  _unsafe_ascomplex_33(A)

N = 3: (rrr, irr, rir, iir, rri, iri, rii, iii)
im3: (rrr, irr, rir, iir, rri, iri, rii, iii) -> order (1,5,2,6,3,7,4,8), permute 2/5/3, 4/6/7
  -> (rrr, rri, irr, iri, rir, rii, iir, iii) -> (complex, 2, 2) [zrr, zir, zri, zii]
"""
function _unsafe_ascomplex_33!(A::AbstractArray{M}) where M<:Multicomplex{T,3} where {T}
  p = convert(Ptr{T}, pointer(A))
  w = unsafe_wrap(Array{T}, p, (8, length(A)))  # 8 x N
  @inbounds for k = 1:length(A)
      # 2<-5, 5<-3, 3<-2
      tmp = w[2,k]
      w[2,k] = w[5,k]
      w[5,k] = w[3,k]
      w[3,k] = tmp
      # 4<-6, 6<-7, 7<-4
      tmp = w[4,k]
      w[4,k] = w[6,k]
      w[6,k] = w[7,k]
      w[7,k] = tmp
  end
  p2 = convert(Ptr{Complex{T}}, pointer(A))
  unsafe_wrap(Array{Complex{T}}, p2, (2, 2, size(A)...))
end


function unsafe_fromcomplex!(A::AbstractArray{M}, unit::Integer) where M<:Multicomplex{T,N,C} where {T,N,C}
    unit > 0 || throw(ArgumentError("unit must be positive"))
    unit <= N || throw(ArgumentError("unit must be less than multicomplex order"))

    if N == 1
        # no action required
    elseif N == 2
        if unit == 1
            # no action required
        else
            # reverse permutation
            _unsafe_ascomplex_22!(A)
        end
    elseif N == 3
        if unit == 1
            # no action required
        elseif unit == 2
            # reverse permutation
            _unsafe_ascomplex_32!(A)
        else
            # reverse 3-cycle permutation
            _unsafe_fromcomplex_33!(A)
        end
    else
        throw(ArgumentError("unsupported multicomplex order"))
    end
    return A
end

"""
  _unsafe_fromcomplex_33(A)

N = 3: (rrr, irr, rir, iir, rri, iri, rii, iii)
im3: (rrr, irr, rir, iir, rri, iri, rii, iii) -> order (1,5,2,6,3,7,4,8), permute 2/5/3, 4/6/7
  -> (rrr, rri, irr, iri, rir, rii, iir, iii) -> (complex, 2, 2) [zrr, zir, zri, zii]
"""
function _unsafe_fromcomplex_33!(A::AbstractArray{M}) where M<:Multicomplex{T,3} where {T}
  p = convert(Ptr{T}, pointer(A))
  w = unsafe_wrap(Array{T}, p, (8, length(A)))  # 8 x N
  @inbounds for k = 1:length(A)
      # 2->5, 5->3, 3->2
      tmp = w[3,k]
      w[3,k] = w[5,k]
      w[5,k] = w[2,k]
      w[2,k] = tmp
      # 4->6, 6->7, 7->4
      tmp = w[7,k]
      w[7,k] = w[6,k]
      w[6,k] = w[4,k]
      w[4,k] = tmp
  end
  p2 = convert(Ptr{Complex{T}}, pointer(A))
  unsafe_wrap(Array{Complex{T}}, p2, (2, 2, size(A)...))
end