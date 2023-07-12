var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = MulticomplexNumbers","category":"page"},{"location":"api/#Module","page":"API","title":"Module","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MulticomplexNumbers","category":"page"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MulticomplexNumber","category":"page"},{"location":"api/#Functions-and-methods","page":"API","title":"Functions and methods","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Multicomplex{N}(value::SVector{C, T}) where {T, N, C}\nMulticomplex(x::Real)\nMulticomplex(x::Real, y::Real)\nMulticomplex(z::Complex)\nMulticomplex(x::Real, y::Real, u::Real, v::Real)\nMulticomplex(a::Complex, b::Complex)\nMulticomplex(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}\ncomponent(m::Multicomplex, k)\nrealest(m::Multicomplex)\nBase.real(m::Multicomplex{T,N,C}) where {T,N,C}\nmatrep(m::Multicomplex{T,N,C}) where {T,N,C}\nascomplex(A::AbstractArray{M}) where M<:Multicomplex{T} where {T}\nBase.:(*)(a::Multicomplex{T,N,C}, b::Multicomplex{T,N,C}) where {T,N,C}\nconj(m::Multicomplex) = Multicomplex(real(m),-imag(m))","category":"page"},{"location":"api/#MulticomplexNumbers.Multicomplex-Union{Tuple{StaticArraysCore.SVector{C, T}}, Tuple{C}, Tuple{N}, Tuple{T}} where {T, N, C}","page":"API","title":"MulticomplexNumbers.Multicomplex","text":"Multicomplex{N}(value)\n\nDefines a multicomplex number ℂ_N, with base field and components defined in value.\n\nInput:     value::SVector{C, T})\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.Multicomplex-Tuple{Real}","page":"API","title":"MulticomplexNumbers.Multicomplex","text":"Multicomplex{N}(value)\n\nDefines a multicomplex number ℂ_N, with base field and components defined in value.\n\nInput:     value::SVector{C, T})\n\n\n\n\n\nBuild a higher order multicomplex number from two multicomplex 'real' and 'imaginary' components\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.Multicomplex-Tuple{Complex}","page":"API","title":"MulticomplexNumbers.Multicomplex","text":"Multicomplex{N}(value)\n\nDefines a multicomplex number ℂ_N, with base field and components defined in value.\n\nInput:     value::SVector{C, T})\n\n\n\n\n\nBuild a higher order multicomplex number from two multicomplex 'real' and 'imaginary' components\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.Multicomplex-NTuple{4, Real}","page":"API","title":"MulticomplexNumbers.Multicomplex","text":"Multicomplex(x,y,u,v) = x + i¹y + i²u + i¹i²v\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.Multicomplex-Union{Tuple{C}, Tuple{N}, Tuple{T}, Tuple{Multicomplex{T, N, C}, Multicomplex{T, N, C}}} where {T, N, C}","page":"API","title":"MulticomplexNumbers.Multicomplex","text":"Build a higher order multicomplex number from two multicomplex 'real' and 'imaginary' components\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.component-Tuple{Multicomplex, Any}","page":"API","title":"MulticomplexNumbers.component","text":"Utility function to extract a real-valued component from a multicomplex number\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.realest-Tuple{Multicomplex}","page":"API","title":"MulticomplexNumbers.realest","text":"Extract the 'most real' component\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.matrep-Union{Tuple{Multicomplex{T, N, C}}, Tuple{C}, Tuple{N}, Tuple{T}} where {T, N, C}","page":"API","title":"MulticomplexNumbers.matrep","text":"matrix representation of a multicomplex number\n\n\n\n\n\n","category":"method"},{"location":"api/#MulticomplexNumbers.ascomplex-Union{Tuple{AbstractArray{M}}, Tuple{M}, Tuple{T}} where {T, M<:(Multicomplex{T})}","page":"API","title":"MulticomplexNumbers.ascomplex","text":"ascomplex(A::AbstractArray{M})\n\nReturns a view of the multicomplex input array A as an array of complex numbers, mapping i⁽¹⁾ -> im\n\nIf A has size (m, n, ...) with multicomplex numbers of dimension N, then the output will have size (2^(N-1), m, n, ...).\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.:*-Union{Tuple{C}, Tuple{N}, Tuple{T}, Tuple{Multicomplex{T, N, C}, Multicomplex{T, N, C}}} where {T, N, C}","page":"API","title":"Base.:*","text":"*(a,b)\n\nMulticomplex multiplication via the matrix representation.\n\n\n\n\n\n","category":"method"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]\nModule = [\"MulticomplexNumbers.jl\"]\nOrder = [:type, :function]","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"CurrentModule = MulticomplexNumbers","category":"page"},{"location":"examples/#Multicomplex-imaginary-units-square-to-1","page":"Examples","title":"Multicomplex imaginary units square to -1","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Multicomplex units im1 through im6 are defined for convenience. Each imaginary unit squares to -1.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> using MulticomplexNumbers\n\njulia> im1\nMulticomplex{Int8, 1, 2}(Int8[0, 1])\n\njulia> im1*im1\nMulticomplex{Int8, 1, 2}(Int8[-1, 0])\n\njulia> im2\nMulticomplex{Int8, 2, 4}(Int8[0, 0, 1, 0])\n\njulia> im2*im2\nMulticomplex{Int8, 2, 4}(Int8[-1, 0, 0, 0])","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"However, the product of two different complex units does not equal -1, but rather represents terms such as i_1i_2==i_2i_1.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> im1*im2\nMulticomplex{Int8, 2, 4}(Int8[0, 0, 0, 1])\n\njulia> im1*im2 == im2*im1\ntrue","category":"page"},{"location":"userguide/#userguide","page":"User guide","title":"User guide","text":"","category":"section"},{"location":"userguide/","page":"User guide","title":"User guide","text":"","category":"page"},{"location":"userguide/","page":"User guide","title":"User guide","text":"CurrentModule = MulticomplexNumbers","category":"page"},{"location":"userguide/","page":"User guide","title":"User guide","text":"For some simple examples, head over to the examples section. For a detailed guide, keep reading.","category":"page"},{"location":"userguide/","page":"User guide","title":"User guide","text":"MulticomplexNumbers.jl is an implementation of multicomplex numbers. These numbers are a generalisation of complex numbers, recursively defined to contain multiple imaginary numbers, i_1, i_2 etc. A new type is defined, MulticomplexNumber, which describes a multicomplex number of any order, and associated arithmetic operations.","category":"page"},{"location":"userguide/","page":"User guide","title":"User guide","text":"The package is loaded as usual:","category":"page"},{"location":"userguide/","page":"User guide","title":"User guide","text":"using MulticomplexNumbers","category":"page"},{"location":"userguide/#Complex-numbers","page":"User guide","title":"Complex numbers","text":"","category":"section"},{"location":"userguide/","page":"User guide","title":"User guide","text":"...","category":"page"},{"location":"userguide/#Bicomplex-numbers","page":"User guide","title":"Bicomplex numbers","text":"","category":"section"},{"location":"userguide/","page":"User guide","title":"User guide","text":"...","category":"page"},{"location":"userguide/#Tricomplex-numbers","page":"User guide","title":"Tricomplex numbers","text":"","category":"section"},{"location":"userguide/","page":"User guide","title":"User guide","text":"...","category":"page"},{"location":"userguide/#Arithmetic-operations","page":"User guide","title":"Arithmetic operations","text":"","category":"section"},{"location":"userguide/#Matrix-representations-and-multicomplex-multiplication","page":"User guide","title":"Matrix representations and multicomplex multiplication","text":"","category":"section"},{"location":"background/#Background","page":"Background","title":"Background","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"","category":"page"},{"location":"background/#Introduction","page":"Background","title":"Introduction","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"MulticomplexNumbers.jl provides a Julia implementation of multicomplex numbers. These numbers are a generalisation of complex numbers, recursively defined to contain multiple imaginary numbers, i_1, i_2 etc. Unlike Clifford algebras, these numbers commute, i.e. i_1i_2=i_2i_1.","category":"page"},{"location":"background/#Algebra-of-the-multicomplex-numbers","page":"Background","title":"Algebra of the multicomplex numbers","text":"","category":"section"},{"location":"background/#Matrix-representations","page":"Background","title":"Matrix representations","text":"","category":"section"},{"location":"background/#Arithmetic","page":"Background","title":"Arithmetic","text":"","category":"section"},{"location":"background/#Other-operations","page":"Background","title":"Other operations","text":"","category":"section"},{"location":"#MulticomplexNumbers.jl","page":"Home","title":"MulticomplexNumbers.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia package for representing multicomplex numbers and performing multicomplex algebra.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Chris Waudby, UCL School of Pharmacy, London (UK).","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you find this package useful, please cite this repository:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Waudby C. A. (2022). MulticomplexNumbers.jl. https://github.com/waudbygroup/MulticomplexNumbers.jl.","category":"page"},{"location":"#Licence","page":"Home","title":"Licence","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MulticomplexNumbers.jl is licenced under the MIT licence; see LICENSE for the full license text.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MulticomplexNumbers.jl is a registered package, and is installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add(\"MulticomplexNumbers\")","category":"page"},{"location":"#Related-packages","page":"Home","title":"Related packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NMRTools.jl: A simple library for handling NMR data in Julia\nDualNumbers.jl: Julia package for representing dual numbers and for performing dual algebra\nForwardDiff.jl: Forward Mode Automatic Differentiation for Julia\nHyperDualNumbers.jl: Julia implementation of HyperDualNumbers\nCliffordAlgebras.jl: A fast and lightweight Julia package for Clifford and geometric algebras\nGrassmann.jl: ⟨Leibniz-Grassmann-Clifford⟩ differential geometric algebra / multivector simplicial complex","category":"page"},{"location":"","page":"Home","title":"Home","text":"<!– ### Acknowledgments","category":"page"},{"location":"","page":"Home","title":"Home","text":"... –>","category":"page"}]
}
