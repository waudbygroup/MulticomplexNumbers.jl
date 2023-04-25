#######
# I/O #
#######
# function show(io::IO, z::Complex)
#     r, i = reim(z)
#     compact = get(io, :compact, false)
#     show(io, r)
#     if signbit(i) && !isnan(i)
#         print(io, compact ? "-" : " - ")
#         if isa(i,Signed) && !isa(i,BigInt) && i == typemin(typeof(i))
#             show(io, -widen(i))
#         else
#             show(io, -i)
#         end
#     else
#         print(io, compact ? "+" : " + ")
#         show(io, i)
#     end
#     if !(isa(i,Integer) && !isa(i,Bool) || isa(i,AbstractFloat) && isfinite(i))
#         print(io, "*")
#     end
#     print(io, "im")
# end
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
# function read(s::IO, ::Type{Complex{T}}) where T<:Real
#     r = read(s,T)
#     i = read(s,T)
#     Complex{T}(r,i)
# end
# function write(s::IO, z::Complex)
#     write(s,real(z),imag(z))
# end