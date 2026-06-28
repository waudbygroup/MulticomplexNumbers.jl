using MulticomplexNumbers
using Zygote
using StaticArrays
using Test

# Central finite-difference gradient over the components of a multicomplex number.
function fdgrad(f, m::Multicomplex{T,N,C}; h = 1e-6) where {T,N,C}
    g = ntuple(C) do k
        mp = Multicomplex{N}(setindex(m.value, m.value[k] + h, k))
        mm = Multicomplex{N}(setindex(m.value, m.value[k] - h, k))
        (f(mp) - f(mm)) / (2h)
    end
    return Multicomplex{N}(SVector{C}(g))
end

# Extract the component vector from whatever Zygote returns as the gradient.
gcomps(g::Multicomplex) = collect(flat(g))
gcomps(g) = collect(g.value)

# Representative multicomplex inputs (chosen away from singularities of log/inv/tan).
samples = [
    Multicomplex(0.4, 0.3),
    Multicomplex(0.4, 0.3, 0.2, 0.1),
    Multicomplex(0.3, 0.2, 0.1, 0.05, 0.15, 0.25, 0.12, 0.08),
]

# Real-valued scalar functions of a single multicomplex argument.
funcs = Dict(
    "identity-realest" => m -> realest(m),
    "sum-abs2"         => m -> abs2(m),
    "abs"              => m -> abs(m + 1),
    "mul"              => m -> realest(m * m),
    "exp"              => m -> realest(exp(m)),
    "log"              => m -> realest(log(m + 1)),
    "sqrt"             => m -> realest(sqrt(m + 1)),
    "inv"              => m -> realest(inv(m + 1)),
    "div"              => m -> realest((m + 1) / (m + 2)),
    "sin"              => m -> realest(sin(m)),
    "cos"              => m -> realest(cos(m)),
    "sinh"             => m -> realest(sinh(m)),
    "cosh"             => m -> realest(cosh(m)),
    "tan"              => m -> realest(tan(m)),
    "tanh"             => m -> realest(tanh(m)),
    "pow-int"          => m -> realest(m^3),
    "pow-literal"      => m -> realest(m * m * m + m^2),
    "conj"             => m -> realest(conj(m) * m),
    "add-scalar"       => m -> realest(2.5 + m - 1.0),
    "scale"            => m -> realest(3.0 * m / 2.0),
    "component"        => m -> component(m, 1) + realest(m),
    "real-imag"        => m -> realest(real(m) + imag(m)),
    "compose"          => m -> realest(exp(sin(m)) + cosh(m * m)),
)

@testset "single-argument gradients" begin
    for (name, f) in funcs
        @testset "$name" begin
            for m in samples
                g = Zygote.gradient(f, m)[1]
                @test g isa Multicomplex
                @test gcomps(g) ≈ gcomps(fdgrad(f, m)) rtol = 1e-4 atol = 1e-5
            end
        end
    end
end

@testset "gradient is a Multicomplex of the input order" begin
    for m in samples
        g = Zygote.gradient(x -> realest(exp(x)), m)[1]
        @test order(g) == order(m)
        @test length(flat(g)) == length(flat(m))
    end
end

@testset "two-argument gradients" begin
    for (a, b) in [
        (Multicomplex(0.4, 0.3), Multicomplex(-0.2, 0.5)),
        (Multicomplex(0.4, 0.3, 0.2, 0.1), Multicomplex(0.1, -0.2, 0.3, 0.05)),
    ]
        h = 1e-6
        ga, gb = Zygote.gradient((x, y) -> realest(x * y + x / (y + 1)), a, b)
        # finite differences
        f1 = x -> realest(x * b + x / (b + 1))
        f2 = y -> realest(a * y + a / (y + 1))
        @test gcomps(ga) ≈ gcomps(fdgrad(f1, a)) rtol = 1e-4 atol = 1e-5
        @test gcomps(gb) ≈ gcomps(fdgrad(f2, b)) rtol = 1e-4 atol = 1e-5
    end
end

@testset "differentiate through real entry points" begin
    # gradient with respect to a real parameter feeding a multicomplex computation
    f = x -> realest(exp(Multicomplex(x, 0.5)))
    gx = Zygote.gradient(f, 0.7)[1]
    fd = (f(0.7 + 1e-6) - f(0.7 - 1e-6)) / 2e-6
    @test gx ≈ fd rtol = 1e-5

    # N=2 real entry
    g = x -> realest(sin(Multicomplex(x, 0.2, 0.1, 0.3)))
    gx2 = Zygote.gradient(g, 0.4)[1]
    fd2 = (g(0.4 + 1e-6) - g(0.4 - 1e-6)) / 2e-6
    @test gx2 ≈ fd2 rtol = 1e-5
end

@testset "differentiate through pair constructor" begin
    # Multicomplex(a, b) builds a higher-order number from two halves
    f = (a, b) -> realest(exp(Multicomplex(a, b)))
    a = Multicomplex(0.3, 0.1)
    b = Multicomplex(0.2, -0.1)
    ga, gb = Zygote.gradient(f, a, b)
    fa = x -> realest(exp(Multicomplex(x, b)))
    fb = y -> realest(exp(Multicomplex(a, y)))
    @test gcomps(ga) ≈ gcomps(fdgrad(fa, a)) rtol = 1e-4 atol = 1e-5
    @test gcomps(gb) ≈ gcomps(fdgrad(fb, b)) rtol = 1e-4 atol = 1e-5
end
