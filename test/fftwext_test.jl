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
    end
end