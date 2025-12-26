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
0 + 1*im1

julia> im1*im1
-1 + 0*im1

julia> im2
(0 + 0*im1) + (1 + 0*im1)*im2

julia> im2*im2
(-1 + 0*im1) + (0 + 0*im1)*im2
```

However, the product of two different complex units does not equal -1, but rather represents terms such as ``i_1i_2==i_2i_1``.
```jldoctest 1
julia> im1*im2
(0 + 0*im1) + (0 + 1*im1)*im2

julia> im1*im2 == im2*im1
true
```
