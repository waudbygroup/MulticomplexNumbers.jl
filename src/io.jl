#######
# I/O #
#######
function Base.show(io::IO, m::Multicomplex{T,0}) where {T}
    r = real(m)
    show(io, real(m))
end
function Base.show(io::IO, m::Multicomplex{T,1}) where {T}
    r, i = reim(m)
    compact = get(io, :compact, false)
    show(io, r)
    if signbit(i) && !isnan(i)
        print(io, compact ? "-" : " - ")
        if isa(i,Signed) && !isa(i,BigInt) && i == typemin(typeof(i))
            show(io, -widen(i))
        else
            show(io, -i)
        end
    else
        print(io, compact ? "+" : " + ")
        show(io, i)
    end
    # if !(isa(i,Integer) && !isa(i,Bool) || isa(i,AbstractFloat) && isfinite(i))
    print(io, "*")
    # end
    print(io, "im1")
end
function Base.show(io::IO, m::Multicomplex{T,N}) where {T,N}
    r, i = reim(m)
    compact = get(io, :compact, false)
    print(io, "(")
    show(io, r)
    print(io, ")")
    print(io, compact ? "+(" : " + (")
    show(io, i)
    print(io, ")*im")
    print(io, order(m))
end
# show(io::IO, z::Complex{Bool}) =
#     print(io, z == im ? "im" : "Complex($(z.re),$(z.im))")

# function show_unquoted(io::IO, z::Complex, ::Int, prec::Int)
#     if operator_precedence(:+) <= prec
#         print(io, "(")
#         show(io, z)
#         print(io, ")")
#     else
#         show(io, z)
#     end
# end

function Base.read(s::IO, ::Type{Multicomplex{T,N}}) where {T,N}
    m = read(s, SVector{2^N,T})
    Multicomplex{N}(m)
end

function Base.write(s::IO, m::Multicomplex)
    write(s, flat(m))
end

"byte order swaps: components are swapped individually"
Base.bswap(m::Multicomplex{T,N}) where {T,N} = Multicomplex{N}(bswap.(m.value))
