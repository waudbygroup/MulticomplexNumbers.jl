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
    elseif N == 4
        if unit == 1
            _unsafe_ascomplex_41!(A)
        elseif unit == 2
            _unsafe_ascomplex_42!(A)
        elseif unit == 3
            _unsafe_ascomplex_43!(A)
        else
            _unsafe_ascomplex_44!(A)
        end
    else
        throw(ArgumentError(
            "unsafe_ascomplex! does not yet support multicomplex order N=$N. " *
            "Currently supported orders: 1, 2, 3, 4. " *
            "Please file an issue at https://github.com/waudbygroup/MulticomplexNumbers.jl/issues " *
            "if you need support for higher orders."
        ))
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


"""
    _unsafe_ascomplex_41!(A)

N = 4: 16 components
im1: Direct reinterpret to (complex, 2, 2, 2) - no permutation needed
"""
function _unsafe_ascomplex_41!(A::AbstractArray{M}) where M<:Multicomplex{T,4} where {T}
    p = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p, (2, 2, 2, size(A)...))
end

"""
    _unsafe_ascomplex_42!(A)

N = 4: 16 components
im2: Swap positions (2,3), (6,7), (10,11), (14,15) within each 16-element block
"""
function _unsafe_ascomplex_42!(A::AbstractArray{M}) where M<:Multicomplex{T,4} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (4, 4*length(A)))  # 4 x (4*N)
    @inbounds for k = 1:(4*length(A))
        w[2,k], w[3,k] = w[3,k], w[2,k]
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, 2, 2, size(A)...))
end

"""
    _unsafe_ascomplex_43!(A)

N = 4: 16 components
im3: Interleave to pair positions differing in bit 2 (i3)
Permutation within each 8-element block: (1,2,3,4,5,6,7,8) -> (1,5,2,6,3,7,4,8)
"""
function _unsafe_ascomplex_43!(A::AbstractArray{M}) where M<:Multicomplex{T,4} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (8, 2*length(A)))  # 8 x (2*N)
    @inbounds for k = 1:(2*length(A))
        # Apply permutation: 2<-5, 5<-3, 3<-2 and 4<-6, 6<-7, 7<-4
        # Same as _unsafe_ascomplex_33!
        tmp = w[2,k]
        w[2,k] = w[5,k]
        w[5,k] = w[3,k]
        w[3,k] = tmp
        tmp = w[4,k]
        w[4,k] = w[6,k]
        w[6,k] = w[7,k]
        w[7,k] = tmp
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, 2, 2, size(A)...))
end

"""
    _unsafe_ascomplex_44!(A)

N = 4: 16 components
im4: Interleave to pair positions differing in bit 3 (i4)
Permutation: (1,...,16) -> (1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16)
"""
function _unsafe_ascomplex_44!(A::AbstractArray{M}) where M<:Multicomplex{T,4} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (16, length(A)))  # 16 x N
    @inbounds for k = 1:length(A)
        # Interleave first 8 with second 8
        # Target: 1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16
        # Use 8 swaps in a specific order to achieve the shuffle
        # This is a complex permutation, we'll use a temporary buffer approach
        t2,t3,t4,t5,t6,t7,t8 = w[2,k],w[3,k],w[4,k],w[5,k],w[6,k],w[7,k],w[8,k]
        t9,t10,t11,t12,t13,t14,t15 = w[9,k],w[10,k],w[11,k],w[12,k],w[13,k],w[14,k],w[15,k]
        # Position mapping: new[i] = old[perm[i]]
        # 1->1, 2->9, 3->2, 4->10, 5->3, 6->11, 7->4, 8->12
        # 9->5, 10->13, 11->6, 12->14, 13->7, 14->15, 15->8, 16->16
        w[2,k] = t9
        w[3,k] = t2
        w[4,k] = t10
        w[5,k] = t3
        w[6,k] = t11
        w[7,k] = t4
        w[8,k] = t12
        w[9,k] = t5
        w[10,k] = t13
        w[11,k] = t6
        w[12,k] = t14
        w[13,k] = t7
        w[14,k] = t15
        w[15,k] = t8
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, 2, 2, size(A)...))
end


"""
    _unsafe_fromcomplex_43!(A)

Reverse permutation of _unsafe_ascomplex_43!
"""
function _unsafe_fromcomplex_43!(A::AbstractArray{M}) where M<:Multicomplex{T,4} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (8, 2*length(A)))  # 8 x (2*N)
    @inbounds for k = 1:(2*length(A))
        # Reverse: 2->5->3->2 and 4->6->7->4
        tmp = w[3,k]
        w[3,k] = w[5,k]
        w[5,k] = w[2,k]
        w[2,k] = tmp
        tmp = w[7,k]
        w[7,k] = w[6,k]
        w[6,k] = w[4,k]
        w[4,k] = tmp
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, 2, 2, size(A)...))
end

"""
    _unsafe_fromcomplex_44!(A)

Reverse permutation of _unsafe_ascomplex_44!
"""
function _unsafe_fromcomplex_44!(A::AbstractArray{M}) where M<:Multicomplex{T,4} where {T}
    p = convert(Ptr{T}, pointer(A))
    w = unsafe_wrap(Array{T}, p, (16, length(A)))  # 16 x N
    @inbounds for k = 1:length(A)
        # Reverse the interleaving
        # Original was: 1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16
        # Need to restore: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
        t2,t3,t4,t5,t6,t7,t8 = w[2,k],w[3,k],w[4,k],w[5,k],w[6,k],w[7,k],w[8,k]
        t9,t10,t11,t12,t13,t14,t15 = w[9,k],w[10,k],w[11,k],w[12,k],w[13,k],w[14,k],w[15,k]
        # Inverse mapping
        w[2,k] = t3
        w[3,k] = t5
        w[4,k] = t7
        w[5,k] = t9
        w[6,k] = t11
        w[7,k] = t13
        w[8,k] = t15
        w[9,k] = t2
        w[10,k] = t4
        w[11,k] = t6
        w[12,k] = t8
        w[13,k] = t10
        w[14,k] = t12
        w[15,k] = t14
    end
    p2 = convert(Ptr{Complex{T}}, pointer(A))
    unsafe_wrap(Array{Complex{T}}, p2, (2, 2, 2, size(A)...))
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
    elseif N == 4
        if unit == 1
            # no action required
        elseif unit == 2
            # reverse permutation (same operation reverses itself)
            _unsafe_ascomplex_42!(A)
        elseif unit == 3
            # reverse 3-cycle permutation
            _unsafe_fromcomplex_43!(A)
        else
            # reverse interleaving
            _unsafe_fromcomplex_44!(A)
        end
    else
        throw(ArgumentError(
            "unsafe_fromcomplex! does not yet support multicomplex order N=$N. " *
            "Currently supported orders: 1, 2, 3, 4. " *
            "Please file an issue at https://github.com/waudbygroup/MulticomplexNumbers.jl/issues " *
            "if you need support for higher orders."
        ))
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