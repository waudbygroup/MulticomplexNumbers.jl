# Executable versions of the examples shown in the README and documentation,
# so that the advertised usage cannot silently break. Kept deliberately small
# and dependency-free (no FFTW); the FFT examples are covered by fftwext_test.jl.

using MulticomplexNumbers
using Test

@testset "README quick start" begin
    z = 1.0 + 2.0 * im1           # order 1 (complex)
    w = 1.0 + 2.0 * im1 + 3.0 * im2   # order 2 (bicomplex)

    @test order(z) == 1
    @test order(w) == 2

    # arithmetic runs and promotes to the higher order
    @test order(z * w) == 2
    @test w + w ≈ 2.0 * w
    @test exp(w) isa Multicomplex
    @test sqrt(w)^2 ≈ w
end

@testset "multicomplex-step differentiation" begin
    f(x) = x^3

    # First derivative: f'(2) = 3*2^2 = 12, recovered with no cancellation.
    h = 1e-100
    x = 2.0 + h * im1
    @test imag(f(x)) / h ≈ 12.0

    # Second derivative via the bicomplex step: f''(2) = 6*2 = 12.
    # The i1*i2 component (index 4) holds h^2 * f''(x).
    y = 2.0 + h * im1 + h * im2
    @test component(f(y), 4) / h^2 ≈ 12.0
end
