using MulticomplexNumbers
using StaticArrays
using Test

@testset "MulticomplexNumbers.jl - constructors" begin
    @inferred Multicomplex(1.)
    m0 = Multicomplex(1.)
    m1 = Multicomplex(m0, m0)
    m2 = Multicomplex(m1, m1)
    m3 = Multicomplex{3}(SVector{8}(1:8))
    m4 = Multicomplex{4}(SVector{16}(1:16))
    
    @inferred Multicomplex(1, 2)
    @inferred Multicomplex(1, 2.5)
    @inferred Multicomplex(true, 2.5)
    @inferred Multicomplex(1 + 2im)
    @inferred Multicomplex{1}(SVector{2}(1:2))
    @inferred Multicomplex(m0, m0)
    @test Multicomplex(1,2) == Multicomplex{1}(SVector{2}(1:2))

    @inferred Multicomplex(1, 2.5, 3, 4)
    @inferred Multicomplex(1+2.5im, 3+4im)
    @inferred Multicomplex(1+2.5im, 3)
    @inferred Multicomplex(1, 3.5+4im)
    @inferred Multicomplex(m2, m2)
    @test Multicomplex(1,2,3,4) == Multicomplex{2}(SVector{4}(1:4))

    @inferred Multicomplex{3}(SVector{8}(1:8))
    @inferred Multicomplex(m2, m2)
    @test Multicomplex(
            Multicomplex{2}(SVector{4}(1:4)),
            Multicomplex{2}(SVector{4}(5:8))
        ) == Multicomplex{3}(SVector{8}(1:8))

    @inferred Multicomplex{4}(SVector{16}(1:16))
    @inferred Multicomplex(m3, m3)
    @test Multicomplex(
            Multicomplex{3}(SVector{8}(1:8)),
            Multicomplex{3}(SVector{8}(9:16))
        ) == Multicomplex{4}(SVector{16}(1:16))

    @test Multicomplex(m0) == m0
    @test Multicomplex(m1) == m1
    @test Multicomplex(m2) == m2
    @test Multicomplex(m3) == m3
    @test Multicomplex(m4) == m4
end

@testset "Real/imag components" begin
    m0 = Multicomplex{0}(SVector{1}(1))
    m1 = Multicomplex{1}(SVector{2}(1:2))
    m2 = Multicomplex{2}(SVector{4}(1:4))
    m3 = Multicomplex{3}(SVector{8}(1:8))
    m4 = Multicomplex{4}(SVector{16}(1:16))

    # order
    @test order(m0) == 0
    @test order(m1) == 1
    @test order(m2) == 2
    @test order(m3) == 3
    @test order(m4) == 4

    # flat
    @test flat(m0) == [1]
    @test flat(m1) == [1,2]
    @test flat(m2) == [1,2,3,4]
    @test flat(m3) == [1,2,3,4,5,6,7,8]
    @test flat(m4) == [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

    #realest
    @test realest(m0) == 1
    @test realest(m1) == 1
    @test realest(m2) == 1
    @test realest(m3) == 1
    @test realest(m4) == 1

    #real
    @test real(m0) == 1
    @test real(m1) == 1
    @test real(m2) == m1
    @test real(m3) == m2
    @test real(m4) == m3

    #imag
    @test imag(m0) == 0
    @test imag(m1) == 2
    @test imag(m2) == Multicomplex{1}(SVector{2}(3:4))
    @test imag(m3) == Multicomplex{2}(SVector{4}(5:8))
    @test imag(m4) == Multicomplex{3}(SVector{8}(9:16))
    
    #reim
    @test reim(m0) == (1, 0)
    @test reim(m1) == (1, 2)
    @test reim(m2) == (m1, Multicomplex{1}(SVector{2}(3:4)))
    @test reim(m3) == (m2, Multicomplex{2}(SVector{4}(5:8)))
    @test reim(m4) == (m3, Multicomplex{3}(SVector{8}(9:16)))

    #isreal
    @test isreal(m0) == true
    @test isreal(m1) == false
    @test isreal(m2) == false
    @test isreal(Multicomplex{1}(SVector{2}(1,0))) == true
    @test isreal(Multicomplex{2}(SVector{4}(1,0,0,0))) == true
    @test isreal(Multicomplex{3}(SVector{8}(1,0,0,0,0,0,0,0))) == true
    
end

@testset "Matrix representations" begin
    m0 = Multicomplex{0}(SVector{1}(1))
    m1 = Multicomplex{1}(SVector{2}(1:2))
    m2 = Multicomplex{2}(SVector{4}(1:4))
    m3 = Multicomplex{3}(SVector{8}(1:8))
    m4 = Multicomplex{4}(SVector{16}(1:16))
    # m5 = Multicomplex{5}(SVector{32}(1:32))

    # matrep
    @test matrep(m0) == SMatrix{1,1}(1)
    @test matrep(m1) == SMatrix{2,2}(1,2,-2,1)
    @test matrep(m2) == SMatrix{4,4}(1, 2, 3, 4, -2, 1, -4, 3, -3, -4, 1, 2, 4, -3, -2, 1)
    @test matrep(m3) == SMatrix{8,8}(1, 2, 3, 4, 5, 6, 7, 8, -2, 1, -4, 3, -6, 5, -8, 7, -3, -4, 1, 2, -7, -8, 5, 6, 4, -3, -2, 1, 8, -7, -6, 5, -5, -6, -7, -8, 1, 2, 3, 4, 6, -5, 8, -7, -2, 1, -4, 3, 7, 8, -5, -6, -3, -4, 1, 2, -8, 7, 6, -5, 4, -3, -2, 1)
    @test matrep(m4) == SMatrix{16,16}(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, -2, 1, -4, 3, -6, 5, -8, 7, -10, 9, -12, 11, -14, 13, -16, 15, -3, -4, 1, 2, -7, -8, 5, 6, -11, -12, 9, 10, -15, -16, 13, 14, 4, -3, -2, 1, 8, -7, -6, 5, 12, -11, -10, 9, 16, -15, -14, 13, -5, -6, -7, -8, 1, 2, 3, 4, -13, -14, -15, -16, 9, 10, 11, 12, 6, -5, 8, -7, -2, 1, -4, 3, 14, -13, 16, -15, -10, 9, -12, 11, 7, 8, -5, -6, -3, -4, 1, 2, 15, 16, -13, -14, -11, -12, 9, 10, -8, 7, 6, -5, 4, -3, -2, 1, -16, 15, 14, -13, 12, -11, -10, 9, -9, -10, -11, -12, -13, -14, -15, -16, 1, 2, 3, 4, 5, 6, 7, 8, 10, -9, 12, -11, 14, -13, 16, -15, -2, 1, -4, 3, -6, 5, -8, 7, 11, 12, -9, -10, 15, 16, -13, -14, -3, -4, 1, 2, -7, -8, 5, 6, -12, 11, 10, -9, -16, 15, 14, -13, 4, -3, -2, 1, 8, -7, -6, 5, 13, 14, 15, 16, -9, -10, -11, -12, -5, -6, -7, -8, 1, 2, 3, 4, -14, 13, -16, 15, 10, -9, 12, -11, 6, -5, 8, -7, -2, 1, -4, 3, -15, -16, 13, 14, 11, 12, -9, -10, 7, 8, -5, -6, -3, -4, 1, 2, 16, -15, -14, 13, -12, 11, 10, -9, -8, 7, 6, -5, 4, -3, -2, 1)
    # @test matrep(m5) == SMatrix{32,32}(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, -2, 1, -4, 3, -6, 5, -8, 7, -10, 9, -12, 11, -14, 13, -16, 15, -18, 17, -20, 19, -22, 21, -24, 23, -26, 25, -28, 27, -30, 29, -32, 31, -3, -4, 1, 2, -7, -8, 5, 6, -11, -12, 9, 10, -15, -16, 13, 14, -19, -20, 17, 18, -23, -24, 21, 22, -27, -28, 25, 26, -31, -32, 29, 30, 4, -3, -2, 1, 8, -7, -6, 5, 12, -11, -10, 9, 16, -15, -14, 13, 20, -19, -18, 17, 24, -23, -22, 21, 28, -27, -26, 25, 32, -31, -30, 29, -5, -6, -7, -8, 1, 2, 3, 4, -13, -14, -15, -16, 9, 10, 11, 12, -21, -22, -23, -24, 17, 18, 19, 20, -29, -30, -31, -32, 25, 26, 27, 28, 6, -5, 8, -7, -2, 1, -4, 3, 14, -13, 16, -15, -10, 9, -12, 11, 22, -21, 24, -23, -18, 17, -20, 19, 30, -29, 32, -31, -26, 25, -28, 27, 7, 8, -5, -6, -3, -4, 1, 2, 15, 16, -13, -14, -11, -12, 9, 10, 23, 24, -21, -22, -19, -20, 17, 18, 31, 32, -29, -30, -27, -28, 25, 26, -8, 7, 6, -5, 4, -3, -2, 1, -16, 15, 14, -13, 12, -11, -10, 9, -24, 23, 22, -21, 20, -19, -18, 17, -32, 31, 30, -29, 28, -27, -26, 25, -9, -10, -11, -12, -13, -14, -15, -16, 1, 2, 3, 4, 5, 6, 7, 8, -25, -26, -27, -28, -29, -30, -31, -32, 17, 18, 19, 20, 21, 22, 23, 24, 10, -9, 12, -11, 14, -13, 16, -15, -2, 1, -4, 3, -6, 5, -8, 7, 26, -25, 28, -27, 30, -29, 32, -31, -18, 17, -20, 19, -22, 21, -24, 23, 11, 12, -9, -10, 15, 16, -13, -14, -3, -4, 1, 2, -7, -8, 5, 6, 27, 28, -25, -26, 31, 32, -29, -30, -19, -20, 17, 18, -23, -24, 21, 22, -12, 11, 10, -9, -16, 15, 14, -13, 4, -3, -2, 1, 8, -7, -6, 5, -28, 27, 26, -25, -32, 31, 30, -29, 20, -19, -18, 17, 24, -23, -22, 21, 13, 14, 15, 16, -9, -10, -11, -12, -5, -6, -7, -8, 1, 2, 3, 4, 29, 30, 31, 32, -25, -26, -27, -28, -21, -22, -23, -24, 17, 18, 19, 20, -14, 13, -16, 15, 10, -9, 12, -11, 6, -5, 8, -7, -2, 1, -4, 3, -30, 29, -32, 31, 26, -25, 28, -27, 22, -21, 24, -23, -18, 17, -20, 19, -15, -16, 13, 14, 11, 12, -9, -10, 7, 8, -5, -6, -3, -4, 1, 2, -31, -32, 29, 30, 27, 28, -25, -26, 23, 24, -21, -22, -19, -20, 17, 18, 16, -15, -14, 13, -12, 11, 10, -9, -8, 7, 6, -5, 4, -3, -2, 1, 32, -31, -30, 29, -28, 27, 26, -25, -24, 23, 22, -21, 20, -19, -18, 17, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30, -31, -32, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, -17, 20, -19, 22, -21, 24, -23, 26, -25, 28, -27, 30, -29, 32, -31, -2, 1, -4, 3, -6, 5, -8, 7, -10, 9, -12, 11, -14, 13, -16, 15, 19, 20, -17, -18, 23, 24, -21, -22, 27, 28, -25, -26, 31, 32, -29, -30, -3, -4, 1, 2, -7, -8, 5, 6, -11, -12, 9, 10, -15, -16, 13, 14, -20, 19, 18, -17, -24, 23, 22, -21, -28, 27, 26, -25, -32, 31, 30, -29, 4, -3, -2, 1, 8, -7, -6, 5, 12, -11, -10, 9, 16, -15, -14, 13, 21, 22, 23, 24, -17, -18, -19, -20, 29, 30, 31, 32, -25, -26, -27, -28, -5, -6, -7, -8, 1, 2, 3, 4, -13, -14, -15, -16, 9, 10, 11, 12, -22, 21, -24, 23, 18, -17, 20, -19, -30, 29, -32, 31, 26, -25, 28, -27, 6, -5, 8, -7, -2, 1, -4, 3, 14, -13, 16, -15, -10, 9, -12, 11, -23, -24, 21, 22, 19, 20, -17, -18, -31, -32, 29, 30, 27, 28, -25, -26, 7, 8, -5, -6, -3, -4, 1, 2, 15, 16, -13, -14, -11, -12, 9, 10, 24, -23, -22, 21, -20, 19, 18, -17, 32, -31, -30, 29, -28, 27, 26, -25, -8, 7, 6, -5, 4, -3, -2, 1, -16, 15, 14, -13, 12, -11, -10, 9, 25, 26, 27, 28, 29, 30, 31, 32, -17, -18, -19, -20, -21, -22, -23, -24, -9, -10, -11, -12, -13, -14, -15, -16, 1, 2, 3, 4, 5, 6, 7, 8, -26, 25, -28, 27, -30, 29, -32, 31, 18, -17, 20, -19, 22, -21, 24, -23, 10, -9, 12, -11, 14, -13, 16, -15, -2, 1, -4, 3, -6, 5, -8, 7, -27, -28, 25, 26, -31, -32, 29, 30, 19, 20, -17, -18, 23, 24, -21, -22, 11, 12, -9, -10, 15, 16, -13, -14, -3, -4, 1, 2, -7, -8, 5, 6, 28, -27, -26, 25, 32, -31, -30, 29, -20, 19, 18, -17, -24, 23, 22, -21, -12, 11, 10, -9, -16, 15, 14, -13, 4, -3, -2, 1, 8, -7, -6, 5, -29, -30, -31, -32, 25, 26, 27, 28, 21, 22, 23, 24, -17, -18, -19, -20, 13, 14, 15, 16, -9, -10, -11, -12, -5, -6, -7, -8, 1, 2, 3, 4, 30, -29, 32, -31, -26, 25, -28, 27, -22, 21, -24, 23, 18, -17, 20, -19, -14, 13, -16, 15, 10, -9, 12, -11, 6, -5, 8, -7, -2, 1, -4, 3, 31, 32, -29, -30, -27, -28, 25, 26, -23, -24, 21, 22, 19, 20, -17, -18, -15, -16, 13, 14, 11, 12, -9, -10, 7, 8, -5, -6, -3, -4, 1, 2, -32, 31, 30, -29, 28, -27, -26, 25, 24, -23, -22, 21, -20, 19, 18, -17, 16, -15, -14, 13, -12, 11, 10, -9, -8, 7, 6, -5, 4, -3, -2, 1)
end

@testset "Arithmetic" begin
    m0 = Multicomplex{0}(SVector{1}(1))
    m1 = Multicomplex{1}(SVector{2}(1:2))
    m2 = Multicomplex{2}(SVector{4}(1:4))
    m3 = Multicomplex{3}(SVector{8}(1:8))
    m4 = Multicomplex{4}(SVector{16}(1:16))
   
    # unary plus and minus
    @test +m2 == Multicomplex{2}(SVector{4}(1,2,3,4))
    @test -m2 == Multicomplex{2}(SVector{4}(-1,-2,-3,-4))

    # addition and subtraction
    @test m0 - m1 == Multicomplex{1}(SVector{2}(0, -2))
    @test m0 + m1 == Multicomplex{1}(SVector{2}(2, 2))
    @test m0 - m1 == Multicomplex{1}(SVector{2}(0, -2))
    @test m0 + m2 == Multicomplex{2}(SVector{4}(2,2,3,4))
    @test m1 + m2 == Multicomplex{2}(SVector{4}(2,4,3,4))

    # multiplication by scalar
    @test 2m0 == Multicomplex{0}(SVector{1}(2))
    @test 2m1 == Multicomplex{1}(SVector{2}(2,4))
    @test 2m2 == Multicomplex{2}(SVector{4}(2,4,6,8))

    # division by scalar
    @test m0/2 == Multicomplex{0}(SVector{1}(0.5))
    @test m1/2 == Multicomplex{1}(SVector{2}(0.5,1))
    @test m2/2 == Multicomplex{2}(SVector{4}(0.5, 1, 1.5, 2))

    # forming combinations
    @test 0.1 + 0.2im1 == Multicomplex{1}(SVector{2}(0.1, 0.2))
    @test 0.1 + 0.2im1 + 0.3im2 + 0.4im1*im2 == Multicomplex{2}(SVector{4}(0.1, 0.2, 0.3, 0.4))
    @test 0.1 + 0.2im1 + 0.3im2 + 0.4im1*im2 + 0.5im3 + 0.6im1*im3 + 0.7im2*im3 + 0.8im1*im2*im3 == Multicomplex{3}(SVector{8}(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))

    # imaginary units square to -1
    @test im1*im1 == Multicomplex(-1)
    @test im2*im2 == Multicomplex(-1)
    @test im3*im3 == Multicomplex(-1)
    @test im4*im4 == Multicomplex(-1)
    @test im5*im5 == Multicomplex(-1)
    @test im6*im6 == Multicomplex(-1)

    # multiplication
    @test m0*m0 == Multicomplex(1)
    @test m1*m1 == Multicomplex(-3+4im)
    @test m0*m1 == Multicomplex(1,2)
    @test m2*m2 == Multicomplex(4-20im, -10+20im)

    # division
    @test (im1/im2)*im2==im1
    @test (im1/im3)*im3==im1
    @test (im2/im3)*im3==im2
    @test (im4/im3)*im3==im4

    # matrix exponential
    @test exp(Multicomplex(1)) == ℯ
    @test exp(0im1) ≈ 1
    @test exp(π*im1) ≈ -1
    @test exp(π*im1/2) ≈ im1
    @test exp(0im2) ≈ 1
    @test exp(π*im2) ≈ -1
    @test exp(π*im2/2) ≈ im2
    @test exp(0im3) ≈ 1
    @test exp(π*im3) ≈ -1
    @test exp(π*im3/2) ≈ im3
end

@testset "Equality" begin
    m0 = Multicomplex{0}(SVector{1}(1))
    m1 = Multicomplex{1}(SVector{2}(1:2))
    m2 = Multicomplex{2}(SVector{4}(1:4))

    @test m0 == m0
    @test m0 == Multicomplex(1)
    @test m1 == m1
    @test m1 == Multicomplex(1+2im)
    @test m2 == m2
    @test m2 == Multicomplex(1+2im, 3+4im)
    @test m0 ≠ m1
    @test m1 ≠ m2

    @test m0 == 1
    @test m1 == 1 + 2im
end


@testset "I/O" begin
    # show
    io = IOBuffer()
    show(io, im1)
    @test String(take!(io)) == "0 + 1*im1"

    io = IOBuffer()
    show(io, 0.5im2)
    @test String(take!(io)) == "(0.0 + 0.0*im1) + (0.5 + 0.0*im1)*im2"

    io = IOBuffer()
    show(io, 1 - 2im1 - 0.5im2 + 0.3im1*im2)
    @test String(take!(io)) == "(1.0 - 2.0*im1) + (-0.5 + 0.3*im1)*im2"

    # int
    m0 = Multicomplex{0}(SVector{1}(1))
    m1 = Multicomplex{1}(SVector{2}(1:2))
    m2 = Multicomplex{2}(SVector{4}(1:4))
    m3 = Multicomplex{3}(SVector{8}(1:8))

    io = IOBuffer()
    write(io, m0)
    write(io, m1)
    write(io, m2)
    write(io, m3)
    seekstart(io)
    m0r = read(io,Multicomplex{Int64,0})
    m1r = read(io,Multicomplex{Int64,1})
    m2r = read(io,Multicomplex{Int64,2})
    m3r = read(io,Multicomplex{Int64,3})

    @test m0r == m0
    @test m1r == m1
    @test m2r == m2
    @test m3r == m3

    # float
    m0 = 0.5 * Multicomplex{0}(SVector{1}(1))
    m1 = 0.5 * Multicomplex{1}(SVector{2}(1:2))
    m2 = 0.5 * Multicomplex{2}(SVector{4}(1:4))
    m3 = 0.5 * Multicomplex{3}(SVector{8}(1:8))

    io = IOBuffer()
    write(io, m0)
    write(io, m1)
    write(io, m2)
    write(io, m3)
    seekstart(io)
    m0r = read(io,Multicomplex{Float64,0})
    m1r = read(io,Multicomplex{Float64,1})
    m2r = read(io,Multicomplex{Float64,2})
    m3r = read(io,Multicomplex{Float64,3})
    
    @test m0r == m0
    @test m1r == m1
    @test m2r == m2
    @test m3r == m3
end

@testset "ascomplex" begin
    @test all(ascomplex([1+2im1]) .== [
        1 + 2im
    ])
    @test all(ascomplex([1+2im1],1) .== [
        1 + 2im
    ])
    @test all(ascomplex([1+2im1+3im2],1) .== [
        1 + 2im
        3 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2],2) .== [
        1 + 3im
        2 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3],1) .== [
        1 + 2im
        3 + 0im
        4 + 0im
        0 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3],2) .== [
        1 + 3im
        2 + 0im
        4 + 0im
        0 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3+5im4]) .== [
        1 + 5im
        2 + 0im
        3 + 0im
        0 + 0im
        4 + 0im
        0 + 0im
        0 + 0im
        0 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3+5im4],1) .== [
        1 + 2im
        3 + 0im
        4 + 0im
        0 + 0im
        5 + 0im
        0 + 0im
        0 + 0im
        0 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3+5im4],2) .== [
        1 + 3im
        2 + 0im
        4 + 0im
        0 + 0im
        5 + 0im
        0 + 0im
        0 + 0im
        0 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3+5im4],3) .== [
        1 + 4im
        2 + 0im
        3 + 0im
        0 + 0im
        5 + 0im
        0 + 0im
        0 + 0im
        0 + 0im
    ])
    @test all(ascomplex([1+2im1+3im2+4im3+5im4],4) .== [
        1 + 5im
        2 + 0im
        3 + 0im
        0 + 0im
        4 + 0im
        0 + 0im
        0 + 0im
        0 + 0im
    ])
end

@testset "zero and one" begin
    # Test zero
    @test @inferred(zero(Multicomplex{Float64,0,1})) == Multicomplex(0.0)
    @test @inferred(zero(Multicomplex{Float64,1,2})) == Multicomplex(0.0, 0.0)
    @test @inferred(zero(Multicomplex{Float64,2,4})) == Multicomplex(0.0, 0.0, 0.0, 0.0)
    @test @inferred(zero(Multicomplex{Int,1,2})) == Multicomplex(0, 0)

    # Test zero from instance
    m2 = Multicomplex(1.0, 2.0, 3.0, 4.0)
    @test @inferred(zero(m2)) == Multicomplex(0.0, 0.0, 0.0, 0.0)

    # Test one
    @test @inferred(one(Multicomplex{Float64,0,1})) == Multicomplex(1.0)
    @test @inferred(one(Multicomplex{Float64,1,2})) == Multicomplex(1.0, 0.0)
    @test @inferred(one(Multicomplex{Float64,2,4})) == Multicomplex(1.0, 0.0, 0.0, 0.0)
    @test @inferred(one(Multicomplex{Int,1,2})) == Multicomplex(1, 0)

    # Test one from instance
    @test @inferred(one(m2)) == Multicomplex(1.0, 0.0, 0.0, 0.0)

    # Test properties: m + zero(m) == m and m * one(m) == m
    m = Multicomplex(1.0, 2.0, 3.0, 4.0)
    @test m + zero(m) == m
    @test m * one(m) == m
    @test one(m) * m == m
end

@testset "fold and isabient" begin
    # Test fold for N=0 (real numbers) - squaring
    m0 = Multicomplex(3.0)
    @test fold(m0) == Multicomplex(9.0)
    @test order(fold(m0)) == 0

    # Test fold for N=1 (complex numbers) - reduces to N=0 (real)
    m1 = Multicomplex(3.0, 4.0)  # 3 + 4i₁
    @test fold(m1) ≈ Multicomplex(25.0)  # |z|² = 3² + 4² = 25
    @test order(fold(m1)) == 0

    m1b = 1.0 + 1.0im1
    @test fold(m1b) ≈ Multicomplex(2.0)  # |1+i|² = 1² + 1² = 2
    @test order(fold(m1b)) == 0

    # Test fold for N=2 - reduces to N=1
    m2 = Multicomplex(1.0, 2.0, 3.0, 4.0)  # 1 + 2i₁ + 3i₂ + 4i₁i₂
    folded_m2 = fold(m2)
    @test order(folded_m2) == 1
    # (1 + 2i₁ + 3i₂ + 4i₁i₂)(1 + 2i₁ - 3i₂ - 4i₁i₂)
    # = (1 + 2i₁)² - (3i₂ + 4i₁i₂)²
    # = (1 + 4i₁ - 4) - (9i₂² + 24i₁i₂² + 16i₁²i₂²)
    # = (-3 + 4i₁) - (-9 - 24 - 16) = (-3 + 4i₁) + 49 = 46 + 4i₁
    @test folded_m2 ≈ Multicomplex(46.0, 4.0)

    # Test the classic abient example: (1 + i₁i₂)
    abient_example = 1.0 + im1*im2
    folded_abient = fold(abient_example)
    @test order(folded_abient) == 1
    @test folded_abient ≈ Multicomplex(0.0, 0.0)  # (1 + i₁i₂)(1 - i₁i₂) = 0

    # Test fold for N=3 - reduces to N=2
    m3 = Multicomplex{3}(SVector{8}(1.0:8.0))
    folded_m3 = fold(m3)
    @test order(folded_m3) == 2
    # The fold should give a specific 4-component result
    @test folded_m3 isa Multicomplex{Float64,2,4}

    # Test isabient for non-abient numbers
    @test isabient(Multicomplex(1.0)) == false
    @test isabient(1.0 + 1.0im1) == false
    @test isabient(Multicomplex(1.0, 2.0, 3.0, 4.0)) == false

    # Test isabient for the classic abient example
    @test isabient(1.0 + im1*im2) == true
    @test isabient(2.0 + 2.0im1*im2) == true

    # Test more abient examples
    # For N=2: a + b*i₁i₂ where a = b is abient (zero divisor property)
    @test isabient(1.0 + im2) == false
    @test isabient(1.0 + im1) == false
    @test isabient(im1*im2) == false  # i₁i₂ alone is not abient (folds to -1, then 1)

    # Test isabient with tolerance
    # This is mathematically abient (a = b form), should fold to ≈ 0
    nearly_abient = Multicomplex(1e-10, 0.0, 0.0, 1e-10)
    @test isabient(nearly_abient, atol=1e-8) == true
    @test isabient(nearly_abient, atol=1e-20) == true  # Should be abient even with tight tolerance

    # Test non-abient number with different tolerances
    not_quite_abient = Multicomplex(1.0, 0.0, 0.0, 2.0)  # a ≠ b, not abient (folds to -3)
    @test isabient(not_quite_abient, atol=1e-8) == false

    # Test type inference
    @test @inferred(fold(Multicomplex(1.0, 2.0))) isa Multicomplex
    @test @inferred(isabient(Multicomplex(1.0, 2.0))) isa Bool
end

@testset "Trigonometric functions" begin
    # Test basic trig functions on real multicomplex (should match real trig)
    m0 = Multicomplex(0.5)
    @test sin(m0) ≈ Multicomplex(sin(0.5))
    @test cos(m0) ≈ Multicomplex(cos(0.5))
    @test tan(m0) ≈ Multicomplex(tan(0.5))

    # Test on N=1 (complex-like)
    m1 = Multicomplex(0.5, 0.3)
    z = 0.5 + 0.3im
    @test sin(m1) ≈ Multicomplex(sin(z))
    @test cos(m1) ≈ Multicomplex(cos(z))
    @test tan(m1) ≈ Multicomplex(tan(z))

    # Test Euler's formula: exp(i*x) = cos(x) + i*sin(x)
    # For multicomplex: exp(x*im1) = cos(x) + im1*sin(x)
    x = 0.5
    m = x * im1  # m = 0.5*im1
    @test exp(m) ≈ Multicomplex(cos(x), sin(x))  # exp(0.5*im1) = cos(0.5) + im1*sin(0.5)

    # Test fundamental identity: sin²(x) + cos²(x) = 1
    m = Multicomplex(0.8, 0.6)
    @test sin(m)^2 + cos(m)^2 ≈ Multicomplex(1.0, 0.0)

    # Test reciprocal functions
    m = Multicomplex(0.5, 0.3)
    @test sec(m) ≈ 1 / cos(m)
    @test csc(m) ≈ 1 / sin(m)
    @test cot(m) ≈ 1 / tan(m)

    # Test type inference
    @test @inferred(sin(Multicomplex(1.0, 2.0))) isa Multicomplex
    @test @inferred(cos(Multicomplex(1.0, 2.0))) isa Multicomplex
    @test @inferred(tan(Multicomplex(1.0, 2.0))) isa Multicomplex
end

@testset "Hyperbolic functions" begin
    # Test basic hyperbolic functions on real multicomplex
    m0 = Multicomplex(0.5)
    @test sinh(m0) ≈ Multicomplex(sinh(0.5))
    @test cosh(m0) ≈ Multicomplex(cosh(0.5))
    @test tanh(m0) ≈ Multicomplex(tanh(0.5))

    # Test on N=1 (complex-like)
    m1 = Multicomplex(0.5, 0.3)
    z = 0.5 + 0.3im
    @test sinh(m1) ≈ Multicomplex(sinh(z))
    @test cosh(m1) ≈ Multicomplex(cosh(z))
    @test tanh(m1) ≈ Multicomplex(tanh(z))

    # Test fundamental identity: cosh²(x) - sinh²(x) = 1
    m = Multicomplex(0.8, 0.6)
    @test cosh(m)^2 - sinh(m)^2 ≈ Multicomplex(1.0, 0.0)

    # Test reciprocal functions
    m = Multicomplex(0.5, 0.3)
    @test sech(m) ≈ 1 / cosh(m)
    @test csch(m) ≈ 1 / sinh(m)
    @test coth(m) ≈ 1 / tanh(m)

    # Test relationship to exp: sinh(x) = (exp(x) - exp(-x))/2
    m = Multicomplex(0.5, 0.3)
    @test sinh(m) ≈ (exp(m) - exp(-m)) / 2

    # Test type inference
    @test @inferred(sinh(Multicomplex(1.0, 2.0))) isa Multicomplex
    @test @inferred(cosh(Multicomplex(1.0, 2.0))) isa Multicomplex
    @test @inferred(tanh(Multicomplex(1.0, 2.0))) isa Multicomplex
end

@testset "Inverse trigonometric functions" begin
    # Test inverse trig functions on real multicomplex
    m0 = Multicomplex(0.5)
    @test asin(m0) ≈ Multicomplex(asin(0.5))
    @test acos(m0) ≈ Multicomplex(acos(0.5))
    @test atan(m0) ≈ Multicomplex(atan(0.5))

    # Test on N=1 (complex-like)
    m1 = Multicomplex(0.5, 0.3)
    z = 0.5 + 0.3im
    @test asin(m1) ≈ Multicomplex(asin(z))
    @test acos(m1) ≈ Multicomplex(acos(z))
    @test atan(m1) ≈ Multicomplex(atan(z))

    # Test inverse relationships
    m = Multicomplex(0.4, 0.3)
    @test sin(asin(m)) ≈ m
    @test cos(acos(m)) ≈ m
    @test tan(atan(m)) ≈ m

    # Test type inference
    @test @inferred(asin(Multicomplex(0.5, 0.3))) isa Multicomplex
    @test @inferred(acos(Multicomplex(0.5, 0.3))) isa Multicomplex
    @test @inferred(atan(Multicomplex(0.5, 0.3))) isa Multicomplex
end

@testset "Inverse hyperbolic functions" begin
    # Test inverse hyperbolic functions on real multicomplex
    m0 = Multicomplex(0.5)
    @test asinh(m0) ≈ Multicomplex(asinh(0.5))
    @test acosh(Multicomplex(1.5)) ≈ Multicomplex(acosh(1.5))
    @test atanh(m0) ≈ Multicomplex(atanh(0.5))

    # Test on N=1 (complex-like)
    m1 = Multicomplex(0.5, 0.3)
    z = 0.5 + 0.3im
    @test asinh(m1) ≈ Multicomplex(asinh(z))
    @test atanh(m1) ≈ Multicomplex(atanh(z))

    # Test inverse relationships
    m = Multicomplex(0.5, 0.3)
    @test sinh(asinh(m)) ≈ m
    @test tanh(atanh(m)) ≈ m

    # Test type inference
    @test @inferred(asinh(Multicomplex(0.5, 0.3))) isa Multicomplex
    @test @inferred(atanh(Multicomplex(0.5, 0.3))) isa Multicomplex
end

@testset "Random number generation" begin
    using Random

    # Test rand for different orders
    @test rand(Multicomplex{Float64,0,1}) isa Multicomplex{Float64,0,1}
    @test rand(Multicomplex{Float64,1,2}) isa Multicomplex{Float64,1,2}
    @test rand(Multicomplex{Float64,2,4}) isa Multicomplex{Float64,2,4}
    @test rand(Multicomplex{Float64,3,8}) isa Multicomplex{Float64,3,8}

    # Test that rand produces values in [0, 1) for Float64
    m = rand(Multicomplex{Float64,2,4})
    @test all(0 ≤ c < 1 for c in flat(m))

    # Test rand with explicit RNG
    rng = MersenneTwister(42)
    m1 = rand(rng, Multicomplex{Float64,1,2})
    @test m1 isa Multicomplex{Float64,1,2}

    # Test reproducibility with seeded RNG
    rng1 = MersenneTwister(123)
    rng2 = MersenneTwister(123)
    @test rand(rng1, Multicomplex{Float64,2,4}) == rand(rng2, Multicomplex{Float64,2,4})

    # Test randn for different orders
    @test randn(Multicomplex{Float64,0,1}) isa Multicomplex{Float64,0,1}
    @test randn(Multicomplex{Float64,1,2}) isa Multicomplex{Float64,1,2}
    @test randn(Multicomplex{Float64,2,4}) isa Multicomplex{Float64,2,4}
    @test randn(Multicomplex{Float64,3,8}) isa Multicomplex{Float64,3,8}

    # Test randn with explicit RNG
    rng = MersenneTwister(42)
    m = randn(rng, Multicomplex{Float64,1,2})
    @test m isa Multicomplex{Float64,1,2}

    # Test that randn produces different values than rand
    rng1 = MersenneTwister(42)
    rng2 = MersenneTwister(42)
    m_rand = rand(rng1, Multicomplex{Float64,1,2})
    m_randn = randn(rng2, Multicomplex{Float64,1,2})
    @test m_rand != m_randn

    # Test reproducibility with seeded RNG for randn
    rng1 = MersenneTwister(456)
    rng2 = MersenneTwister(456)
    @test randn(rng1, Multicomplex{Float64,2,4}) == randn(rng2, Multicomplex{Float64,2,4})

    # Test that rand works with integer types
    m_int = rand(Multicomplex{Int,1,2})
    @test m_int isa Multicomplex{Int,1,2}

    # Test type inference
    @test @inferred(rand(Multicomplex{Float64,1,2})) isa Multicomplex
    @test @inferred(randn(Multicomplex{Float64,1,2})) isa Multicomplex
end