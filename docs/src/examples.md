# [Examples](@id Examples)
---

```@meta
CurrentModule = MulticomplexNumbers
```

## Multicomplex imaginary units square to -1
Multicomplex units `im1` through `im6` are defined for convenience. Each imaginary unit squares to -1.
```jldoctest 1
julia> using MulticomplexNumbers

julia> im1
Multicomplex{Int8, 1, 2}(Int8[0, 1])

julia> im1*im1
Multicomplex{Int8, 1, 2}(Int8[-1, 0])

julia> im2
Multicomplex{Int8, 2, 4}(Int8[0, 0, 1, 0])

julia> im2*im2
Multicomplex{Int8, 2, 4}(Int8[-1, 0, 0, 0])
```

However, the product of two different complex units does not equal -1, but rather represents terms such as ``i_1i_2==i_2i_1``.
```jldoctest 1
julia> im1*im2
Multicomplex{Int8, 2, 4}(Int8[0, 0, 0, 1])

julia> im1*im2 == im2*im1
true
```
