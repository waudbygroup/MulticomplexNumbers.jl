using MulticomplexNumbers
using FFTW
using Test

@testset "FFTW" begin
    testRe = rand(512,256)
    testIm = rand(512,256)
    testC = testRe + 1im*testIm
    ansC = fft(testC)
    ansRe = real(ansC)

    test1 = testRe + im1*testIm
    fft!(test1, 1)
    @test ansC ≈ test1

    test21 = testRe + im1*testIm .+ 0im2
    fft!(test21, 1)
    @test ansRe ≈ realest.(test21)

    test22 = testRe + im2*testIm
    fft!(test22, 2)
    @test ansRe ≈ realest.(test22)

    test31 = testRe + im1*testIm .+ 0im3
    fft!(test31, 1)
    @test ansRe ≈ realest.(test31)

    test32 = testRe + im2*testIm .+ 0im3
    fft!(test32, 2)
    @test ansRe ≈ realest.(test32)

    test33 = testRe + im3*testIm
    fft!(test33, 3)
    @test ansRe ≈ realest.(test33)

    # N=4 tests
    test41 = testRe + im1*testIm .+ 0im4
    fft!(test41, 1)
    @test ansRe ≈ realest.(test41)

    test42 = testRe + im2*testIm .+ 0im4
    fft!(test42, 2)
    @test ansRe ≈ realest.(test42)

    test43 = testRe + im3*testIm .+ 0im4
    fft!(test43, 3)
    @test ansRe ≈ realest.(test43)

    test44 = testRe + im4*testIm
    fft!(test44, 4)
    @test ansRe ≈ realest.(test44)

    # selecting a dimension

    for d=1:3
        testRe = rand(128,256,16)
        testIm = rand(128,256,16)
        testC = testRe + 1im*testIm
        ansC = fft(testC, d)
        ansRe = real(ansC)

        test1 = testRe + im1*testIm
        fft!(test1, 1, d)
        @test ansC ≈ test1

        test21 = testRe + im1*testIm .+ 0im2
        fft!(test21, 1, d)
        @test ansRe ≈ realest.(test21)

        test22 = testRe + im2*testIm
        fft!(test22, 2, d)
        @test ansRe ≈ realest.(test22)

        test31 = testRe + im1*testIm .+ 0im3
        fft!(test31, 1, d)
        @test ansRe ≈ realest.(test31)

        test32 = testRe + im2*testIm .+ 0im3
        fft!(test32, 2, d)
        @test ansRe ≈ realest.(test32)

        test33 = testRe + im3*testIm
        fft!(test33, 3, d)
        @test ansRe ≈ realest.(test33)

        # N=4 tests
        test41 = testRe + im1*testIm .+ 0im4
        fft!(test41, 1, d)
        @test ansRe ≈ realest.(test41)

        test42 = testRe + im2*testIm .+ 0im4
        fft!(test42, 2, d)
        @test ansRe ≈ realest.(test42)

        test43 = testRe + im3*testIm .+ 0im4
        fft!(test43, 3, d)
        @test ansRe ≈ realest.(test43)

        test44 = testRe + im4*testIm
        fft!(test44, 4, d)
        @test ansRe ≈ realest.(test44)
    end
end

@testset "ifft! round-trip" begin
    for N in 1:4
        testRe = rand(64, 32)
        testIm = rand(64, 32)
        unit_mc = imN(N)
        A = testRe + unit_mc * testIm
        A_orig = copy(A)

        fft!(A, N)
        ifft!(A, N)
        @test A ≈ A_orig

        # with specific dim
        A2 = copy(A_orig)
        fft!(A2, N, 1)
        ifft!(A2, N, 1)
        @test A2 ≈ A_orig
    end
end

@testset "allocating fft / ifft" begin
    testRe = rand(64, 32)
    testIm = rand(64, 32)

    for N in 1:3
        unit_mc = imN(N)
        A = testRe + unit_mc * testIm
        A_orig = copy(A)

        # fft should not mutate input
        B = fft(A, N)
        @test A == A_orig
        @test B ≈ fft!(copy(A), N)

        # ifft should not mutate input
        C = ifft(B, N)
        @test C ≈ A_orig
    end
end

@testset "bfft!" begin
    testRe = rand(64, 32)
    testIm = rand(64, 32)
    n = length(testRe)

    for N in 1:3
        A = testRe + imN(N) * testIm
        A_orig = copy(A)

        fft!(A, N)
        bfft!(A, N)
        @test A ≈ A_orig .* n

        # allocating bfft
        A2 = fft(A_orig, N)
        B = bfft(A2, N)
        @test B ≈ A_orig .* n
    end
end

@testset "Multicomplex unit dispatch" begin
    testRe = rand(64, 32)
    testIm = rand(64, 32)

    A1 = testRe + im2 * testIm
    A2 = copy(A1)

    fft!(A1, 2)
    fft!(A2, im2)
    @test A1 ≈ A2

    # allocating variant
    A3 = testRe + im3 * testIm
    @test fft(A3, 3) ≈ fft(A3, im3)
end

@testset "fftshift on multicomplex arrays" begin
    A = [Multicomplex(Float64(i), Float64(i+1)) for i in 1:8]
    shifted = fftshift(A)
    @test length(shifted) == length(A)
    @test shifted[1] == A[5]
    @test shifted[5] == A[1]
end
