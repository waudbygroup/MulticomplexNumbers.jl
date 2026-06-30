---
title: 'MulticomplexNumbers.jl: Commutative hypercomplex algebra in Julia for multidimensional NMR and numerical differentiation'
tags:
  - Julia
  - multicomplex numbers
  - hypercomplex algebra
  - numerical differentiation
  - automatic differentiation
  - NMR spectroscopy
authors:
  - name: Christopher A. Waudby
    orcid: 0000-0001-7810-3753
    affiliation: 1
affiliations:
  - name: UCL School of Pharmacy, University College London, London, United Kingdom
    index: 1
date: 30 June 2026
bibliography: paper.bib
---

# Summary

`MulticomplexNumbers.jl` is a Julia package for representing multicomplex numbers
and performing multicomplex algebra. Multicomplex numbers are a commutative
generalisation of the complex numbers, defined recursively to contain an
arbitrary number of imaginary units $i_1, i_2, \ldots, i_n$, each squaring to
$-1$. Unlike quaternions or other Clifford (geometric) algebras, the imaginary
units commute ($i_j i_k = i_k i_j$), which makes the resulting algebra a natural
setting for analytic computation [@Bell2021; @Casado2020].

The package provides a single parametric `Multicomplex{T,N,C}` number type that
integrates fully with Julia's `Number` interface, including promotion and
conversion rules, so that multicomplex numbers can be used as a drop-in scalar
type in generic numerical code. Here, `T` is the base type, `N` the number of
imaginary units, and `C=2^N`. The package implements the full set of elementary
arithmetic and transcendental operations ($+$, $-$, $\times$, $/$, powers,
`exp`, `log`, `sqrt`, and the trigonometric and hyperbolic functions and their
inverses), together with matrix representations, conjugation, fold/division
algorithm, and views that reinterpret arrays of multicomplex numbers as nested
complex arrays. The type is built on stack-allocated `StaticArrays` storage and
is type-stable, so arithmetic is allocation-free and competitive with hand-written
complex code. An optional package extension adds in-place multicomplex fast
Fourier transforms when `FFTW` is loaded.

# Statement of need

The motivating application for this package is **multi-dimensional nuclear
magnetic resonance (NMR) spectroscopy**. Each independently sampled (indirect)
time dimension of an NMR experiment is acquired in quadrature, carrying its own
real/imaginary (cosine/sine) phase pair, so an $n$-dimensional experiment is
naturally and exactly represented by an array of order-$n$ multicomplex numbers,
with the imaginary unit $i_k$ attached to dimension $k$ [@Delsuc1988]. Holding
the data in this representation lets an entire processing pipeline — apodization,
the per-dimension Fourier transforms, and phase correction applied independently
in each dimension — be written as ordinary multicomplex arithmetic, instead of
manually bookkeeping the $2^n$ interleaved real arrays of a "hypercomplex"
dataset. `MulticomplexNumbers.jl` provides this representation together with an
`FFTW`-backed in-place multicomplex FFT, and underpins the NMR data-handling
package `NMRTools.jl`.

The same algebra also enables **high-order numerical differentiation**. The
complex-step derivative approximation,
$f'(x) \approx \operatorname{Im}\,f(x + i h)/h$, computes first derivatives to
machine precision because it avoids the subtractive cancellation that limits
finite-difference schemes, allowing step sizes as small as $h = 10^{-100}$
[@Squire1998; @Lai2008]. Multicomplex numbers extend this idea to arbitrary
order: by promoting a function argument to a multicomplex number with several
independent infinitesimal imaginary perturbations, mixed and higher-order partial
derivatives can be recovered, again without cancellation error, from the
appropriate components of the result [@Lantoine2012; @Bell2021]. This makes the
multicomplex step method attractive whenever exact, high-order derivatives of an
existing scalar function are required and the source can simply be re-run with a
different number type — for example in sensitivity analysis, optimisation, and
the differentiation of thermodynamic functions [@Bell2021].

Multicomplex arithmetic has been implemented before in other languages — a C++
library from NIST [@Bell2021], the MATLAB class of Algorithm 1008 [@Casado2020],
and the C++/Python library `MultiZ` [@AguirreMesa2020] — but none is usable from
Julia, and all fix the base type to double precision. These libraries are applied
chiefly to exact high-order derivatives and sensitivity analysis in engineering,
for example in solid mechanics, fracture, and design optimisation
[@AguirreMesa2020], while the order-two (bicomplex, or *tessarine*) case also
underlies bicomplex function theory and arises in signal processing and bicomplex
quantum mechanics [@LunaElizarraras2015].

Within the Julia ecosystem, the established tools for derivative computation are
based on dual numbers: `ForwardDiff.jl` [@Revels2016] implements forward-mode
automatic differentiation, `DualNumbers.jl` provides first-order dual numbers,
and `HyperDualNumbers.jl` provides hyper-dual numbers for exact second
derivatives. These carry truncated Taylor coefficients rather than forming a
closed algebra, and are specialised for differentiation rather than offering a
general hypercomplex number type. `CliffordAlgebras.jl` implements Clifford and
geometric algebras, but there the generating units anticommute
($e_j e_k = -e_k e_j$); multicomplex units instead commute, giving a different
algebra better suited to data with several independent, factorised phase
dimensions. `MulticomplexNumbers.jl` is, to our knowledge, the first native Julia
implementation of multicomplex algebra, and it differs from this prior art in
three ways that matter for research use: (i) it is generic in the base type, so
the same code runs in `Float64`, `BigFloat`, or any other `Real`; (ii) it is a
first-class `Number` subtype with promotion and conversion rules, so it composes
with Julia's existing array, linear-algebra, and FFT code without modification;
and (iii) it provides a multicomplex FFT for the multi-dimensional NMR
application, with no counterpart in the reference implementations. The package is
used in production by the NMR data-analysis package `NMRTools.jl`.

The package ships with a test suite that checks both numerical
accuracy and type stability across orders, a benchmark suite run in continuous
integration, and documentation comprising a mathematical background, a user
guide, and worked application examples.


# References
