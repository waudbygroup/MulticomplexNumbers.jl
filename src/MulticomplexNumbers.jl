module MulticomplexNumbers

using LinearAlgebra
using StaticArrays
using PrecompileTools

export Multicomplex
export im1, im2, im3, im4, im5, im6
export order
export flat
export component
export realest
export matrep
export ascomplex
export fold, isabient

include("base.jl")
include("io.jl")
include("representations.jl")
include("arithmetic.jl")

include("precompile.jl")

end
