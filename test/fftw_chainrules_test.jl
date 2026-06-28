using MulticomplexNumbers
using FFTW
using Zygote
using StaticArrays
using Test

gcomps(g::Multicomplex) = collect(flat(g))
gcomps(g) = collect(g.value)

# Central finite-difference gradient over every component of every array element.
function fdgrad_arr(f, A; h = 1e-6)
    G = similar(A)
    for idx in eachindex(A)
        m = A[idx]
        N = order(m)
        C = length(flat(m))
        g = ntuple(C) do k
            Ap = copy(A)
            Am = copy(A)
            Ap[idx] = Multicomplex{N}(setindex(m.value, m.value[k] + h, k))
            Am[idx] = Multicomplex{N}(setindex(m.value, m.value[k] - h, k))
            (f(Ap) - f(Am)) / (2h)
        end
        G[idx] = Multicomplex{N}(SVector{C}(g))
    end
    return G
end

function check(f, A)
    G = Zygote.gradient(f, A)[1]
    Gfd = fdgrad_arr(f, A)
    @test G !== nothing
    @test all(gcomps(G[i]) ≈ gcomps(Gfd[i]) for i in eachindex(A))
end

@testset "fft / ifft / bfft gradients (N=1)" begin
    A = [Multicomplex(0.1, 0.2), Multicomplex(0.3, -0.1),
         Multicomplex(-0.2, 0.4), Multicomplex(0.05, 0.15)]

    check(X -> realest(sum(fft(X, 1))), A)
    check(X -> realest(sum(ifft(X, 1))), A)
    check(X -> realest(sum(bfft(X, 1))), A)
    # nonlinear scalar readout to exercise the full chain
    check(X -> abs2(sum(fft(X, 1))), A)
end

@testset "fft gradients along different units (N=2)" begin
    A = [Multicomplex(0.1, 0.2, 0.05, -0.1) for _ in 1:4]
    A = A .+ [Multicomplex(0.01 * k, 0.02 * k, -0.01 * k, 0.03 * k) for k in 1:4]

    check(X -> realest(sum(fft(X, 1))), A)
    check(X -> realest(sum(fft(X, 2))), A)
    check(X -> realest(sum(ifft(X, 2))), A)
end

@testset "fft with explicit dims (2D array)" begin
    A = [Multicomplex(0.1 * i + 0.2 * j, 0.05 * i - 0.1 * j) for i in 1:3, j in 1:2]
    check(X -> realest(sum(fft(X, 1, 1:2))), A)
    check(X -> realest(sum(fft(X, 1, 1))), A)
end

@testset "round-trip differentiability" begin
    A = [Multicomplex(0.1, 0.2), Multicomplex(0.3, -0.1),
         Multicomplex(-0.2, 0.4), Multicomplex(0.05, 0.15)]
    # ifft(fft(A)) == A, so d/dA realest(sum(ifft(fft(A)))) == d/dA realest(sum(A))
    f = X -> realest(sum(ifft(fft(X, 1), 1)))
    G = Zygote.gradient(f, A)[1]
    Gfd = fdgrad_arr(f, A)
    @test all(gcomps(G[i]) ≈ gcomps(Gfd[i]) for i in eachindex(A))
end
