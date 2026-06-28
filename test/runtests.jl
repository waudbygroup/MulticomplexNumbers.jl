using MulticomplexNumbers
using Test
using SafeTestsets

@safetestset "MulticomplexNumbers" begin
    include("base_test.jl")
end

@safetestset "FFTWExt" begin
    include("fftwext_test.jl")
end

@safetestset "ChainRulesCoreExt" begin
    include("chainrules_test.jl")
end

@safetestset "FFTWChainRulesExt" begin
    include("fftw_chainrules_test.jl")
end
