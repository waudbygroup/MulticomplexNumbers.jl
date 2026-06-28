# Automatic Differentiation

MulticomplexNumbers.jl ships [ChainRulesCore](https://github.com/JuliaDiff/ChainRulesCore.jl)
rules as a package extension, so reverse-mode automatic-differentiation frameworks
such as [Zygote](https://github.com/FluxML/Zygote.jl) can differentiate through
multicomplex algebra. The rules load automatically as soon as `ChainRulesCore`
(or a package that depends on it, like Zygote) is available — there is nothing to
import.

```julia
using MulticomplexNumbers, Zygote

f(m) = realest(exp(m) + m * m)

m = Multicomplex(0.4, 0.3, 0.2, 0.1)   # a ℂ₂ number
g = Zygote.gradient(f, m)[1]
```

## What gets differentiated

A `Multicomplex{T,N,C}` is treated as a real vector space of its `C = 2^N`
components. The (co)tangent of a multicomplex number is therefore **another
multicomplex number of the same order**, holding the partial derivatives with
respect to each component. Gradients of real-valued losses with respect to
multicomplex inputs come back as `Multicomplex` values.

Rules are provided for the core algebra:

  * constructors and component access (`real`, `imag`, `conj`, `flat`,
    `component`, `realest`),
  * arithmetic (`+`, `-`, `*`, `/`, `inv`, `^`, scalar scaling),
  * absolute values (`abs`, `abs2`), and
  * the elementary transcendental functions `exp`, `log`, `sqrt`, `sin`, `cos`,
    `tan`, `sinh`, `cosh`, `tanh`.

When both `FFTW` and `ChainRulesCore` are loaded, the allocating FFT operations
(`fft`, `ifft`, `bfft`) become differentiable as well. The in-place variants
(`fft!`, `ifft!`, `bfft!`) mutate their argument and are not supported by
reverse-mode AD.

```julia
using MulticomplexNumbers, FFTW, Zygote

A = [Multicomplex(0.1, 0.2), Multicomplex(0.3, -0.1),
     Multicomplex(-0.2, 0.4), Multicomplex(0.05, 0.15)]

Zygote.gradient(X -> abs2(sum(fft(X, 1))), A)
```

## The adjoint of multicomplex multiplication

For readers interested in the mathematics: for a real-valued loss the
reverse-mode adjoint of "multiply by `y`" is the transpose of `y`'s
multiplication matrix, `matrep(y)ᵀ`. It can be shown that

```math
\operatorname{matrep}(m)^{\mathsf{T}} = \operatorname{matrep}(\tau(m)),
```

where ``\tau`` flips the sign of *every* imaginary unit (it negates each
component whose index contains an odd number of imaginary factors). Because
``\tau`` is a ring automorphism, the pullback of any analytic function ``f`` is

```math
\bar m = \tau\!\big(f'(m)\big)\, \bar z, \qquad z = f(m).
```

For ``N = 1`` this reduces to ordinary complex conjugation, recovering the
familiar complex multiplication adjoint ``\bar x = \bar z\,\overline{y}``.
