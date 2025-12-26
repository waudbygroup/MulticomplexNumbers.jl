# API

---

```@meta
CurrentModule = MulticomplexNumbers
```

## Module
```@docs
MulticomplexNumbers
```

## Types

```@docs
Multicomplex
```

## Functions and methods

```@docs
Multicomplex{N}(value::SVector{C, T}) where {T, N, C}
Multicomplex(x::Real)
Multicomplex(x::Real, y::Real)
Multicomplex(z::Complex)
Multicomplex(x::Real, y::Real, u::Real, v::Real)
Multicomplex(a::Complex, b::Complex)
Multicomplex(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
component(m::Multicomplex, k)
realest(m::Multicomplex)
Base.real(m::Multicomplex{T,N,C}) where {T,N,C}
matrep(m::Multicomplex{T,N,C}) where {T,N,C}
ascomplex(A::AbstractArray{M}) where M<:Multicomplex{T} where {T}
Base.:(*)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}
conj(m::Multicomplex) = Multicomplex(real(m),-imag(m))
```

## Index

```@index
Pages = ["api.md"]
Module = ["MulticomplexNumbers.jl"]
Order = [:type, :function]
```

