# MulticomplexNumbers

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://waudbygroup.github.io/MulticomplexNumbers.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://waudbygroup.github.io/MulticomplexNumbers.jl/dev)
[![CI](https://github.com/waudbygroup/MulticomplexNumbers.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/waudbygroup/MulticomplexNumbers.jl/actions/workflows/Runtests.yml)
[![codecov](https://codecov.io/gh/waudbygroup/MulticomplexNumbers.jl/branch/main/graph/badge.svg?token=V9ND8Y3R8A)](https://codecov.io/gh/waudbygroup/MulticomplexNumbers.jl)

A Julia package for representing multicomplex numbers and performing multicomplex algebra.

## Multicomplex numbers

Multicomplex numbers are a generalisation of complex numbers, recursively defined to contain multiple imaginary numbers, $i_1$, $i_2$ etc. Unlike Clifford algebras, these numbers commute, i.e. $i_1i_2=i_2i_1$.

* NIST report describing multicomplex algebra and computer implementation (for numberical differentiation): https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf
* NIST C++ implementation: https://github.com/usnistgov/multicomplex/blob/master/multicomplex/src/main.cpp
* Casado JMV, Hewson R. Algorithm 1008: Multicomplex Number Class for Matlab, with a Focus on the Accurate Calculation of Small Imaginary Terms for Multicomplex Step Sensitivity Calculations. ACM Trans Math Softw. 2020;46: 1–26. http://dx.doi.org/10.1145/3378542
* Simple C++ implementation: http://tamivox.org/eugene/multicomplex/index.html


## TODO

* FFT
* The comparison `exp(0im1) == 1.` throws a TypeError


## useful links (notes to self for development)

* Creating a package: https://jaantollander.com/post/how-to-create-software-packages-with-julia-language/#package-structure
* Documentation: https://juliadocs.github.io/Documenter.jl/stable/
* Creating a package: https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/11-developing-julia-packages - this has notes on codecov, managing dependencies for the main package and testing, and using TagBot
* Registering package releases with Registrator/TagBot - https://github.com/JuliaComputing/Registrator.jl
* CodeCov Julia example - https://github.com/codecov/example-julia

## notes for coding

`@boundscheck` - fence off code requiring bound checking - https://docs.julialang.org/en/v1/devdocs/boundscheck/
