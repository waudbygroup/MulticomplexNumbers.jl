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