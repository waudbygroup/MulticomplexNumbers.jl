using BenchmarkTools
using MulticomplexNumbers
using StaticArrays

# Create a BenchmarkGroup to hold all benchmarks
const SUITE = BenchmarkGroup()

# Setup some test data
const m0 = Multicomplex(1.5)
const m1 = Multicomplex(1.5, 2.3)
const m2 = Multicomplex(1.5, 2.3, 3.1, 4.2)
const m3 = Multicomplex{3}(SVector{8}(1:8))
const m4 = Multicomplex{4}(SVector{16}(1:16))

const m1_complex = 1.5 + 2.3im  # For comparison

# Constructors
SUITE["constructors"] = BenchmarkGroup()
SUITE["constructors"]["N=0"] = @benchmarkable Multicomplex(1.5)
SUITE["constructors"]["N=1"] = @benchmarkable Multicomplex(1.5, 2.3)
SUITE["constructors"]["N=2"] = @benchmarkable Multicomplex(1.5, 2.3, 3.1, 4.2)
SUITE["constructors"]["N=3"] = @benchmarkable Multicomplex{3}(SVector{8}(1:8))
SUITE["constructors"]["from complex"] = @benchmarkable Multicomplex($m1_complex)

# Basic arithmetic
SUITE["arithmetic"] = BenchmarkGroup()
SUITE["arithmetic"]["addition N=1"] = @benchmarkable $m1 + $m1
SUITE["arithmetic"]["addition N=2"] = @benchmarkable $m2 + $m2
SUITE["arithmetic"]["addition N=3"] = @benchmarkable $m3 + $m3

SUITE["arithmetic"]["subtraction N=1"] = @benchmarkable $m1 - $m1
SUITE["arithmetic"]["subtraction N=2"] = @benchmarkable $m2 - $m2
SUITE["arithmetic"]["subtraction N=3"] = @benchmarkable $m3 - $m3

SUITE["arithmetic"]["scalar multiplication N=1"] = @benchmarkable 2.5 * $m1
SUITE["arithmetic"]["scalar multiplication N=2"] = @benchmarkable 2.5 * $m2
SUITE["arithmetic"]["scalar multiplication N=3"] = @benchmarkable 2.5 * $m3

SUITE["arithmetic"]["scalar division N=1"] = @benchmarkable $m1 / 2.5
SUITE["arithmetic"]["scalar division N=2"] = @benchmarkable $m2 / 2.5
SUITE["arithmetic"]["scalar division N=3"] = @benchmarkable $m3 / 2.5

# Multicomplex multiplication (expensive)
SUITE["multiplication"] = BenchmarkGroup()
SUITE["multiplication"]["N=0"] = @benchmarkable $m0 * $m0
SUITE["multiplication"]["N=1"] = @benchmarkable $m1 * $m1
SUITE["multiplication"]["N=2"] = @benchmarkable $m2 * $m2
SUITE["multiplication"]["N=3"] = @benchmarkable $m3 * $m3
SUITE["multiplication"]["N=4"] = @benchmarkable $m4 * $m4

# Compare with Complex multiplication
SUITE["multiplication"]["Complex (for comparison)"] = @benchmarkable $m1_complex * $m1_complex

# Division
SUITE["division"] = BenchmarkGroup()
SUITE["division"]["N=1"] = @benchmarkable $m1 / $m1
SUITE["division"]["N=2"] = @benchmarkable $m2 / $m2
SUITE["division"]["N=3"] = @benchmarkable $m3 / $m3

# Powers
SUITE["powers"] = BenchmarkGroup()
SUITE["powers"]["square N=1"] = @benchmarkable $m1 ^ 2
SUITE["powers"]["square N=2"] = @benchmarkable $m2 ^ 2
SUITE["powers"]["cube N=1"] = @benchmarkable $m1 ^ 3
SUITE["powers"]["cube N=2"] = @benchmarkable $m2 ^ 3
SUITE["powers"]["real power N=1"] = @benchmarkable $m1 ^ 2.5
SUITE["powers"]["real power N=2"] = @benchmarkable $m2 ^ 2.5

# Exponentials and logarithms
SUITE["transcendental"] = BenchmarkGroup()
SUITE["transcendental"]["exp N=1"] = @benchmarkable exp($m1)
SUITE["transcendental"]["exp N=2"] = @benchmarkable exp($m2)
SUITE["transcendental"]["exp N=3"] = @benchmarkable exp($m3)

SUITE["transcendental"]["log N=1"] = @benchmarkable log($m1)
SUITE["transcendental"]["log N=2"] = @benchmarkable log($m2)
SUITE["transcendental"]["log N=3"] = @benchmarkable log($m3)

SUITE["transcendental"]["sqrt N=1"] = @benchmarkable sqrt($m1)
SUITE["transcendental"]["sqrt N=2"] = @benchmarkable sqrt($m2)
SUITE["transcendental"]["sqrt N=3"] = @benchmarkable sqrt($m3)

# Trigonometric functions
SUITE["trigonometric"] = BenchmarkGroup()
SUITE["trigonometric"]["sin N=1"] = @benchmarkable sin($m1)
SUITE["trigonometric"]["sin N=2"] = @benchmarkable sin($m2)
SUITE["trigonometric"]["cos N=1"] = @benchmarkable cos($m1)
SUITE["trigonometric"]["cos N=2"] = @benchmarkable cos($m2)
SUITE["trigonometric"]["tan N=1"] = @benchmarkable tan($m1)
SUITE["trigonometric"]["tan N=2"] = @benchmarkable tan($m2)

# Hyperbolic functions
SUITE["hyperbolic"] = BenchmarkGroup()
SUITE["hyperbolic"]["sinh N=1"] = @benchmarkable sinh($m1)
SUITE["hyperbolic"]["sinh N=2"] = @benchmarkable sinh($m2)
SUITE["hyperbolic"]["cosh N=1"] = @benchmarkable cosh($m1)
SUITE["hyperbolic"]["cosh N=2"] = @benchmarkable cosh($m2)
SUITE["hyperbolic"]["tanh N=1"] = @benchmarkable tanh($m1)
SUITE["hyperbolic"]["tanh N=2"] = @benchmarkable tanh($m2)

# Inverse trigonometric functions
SUITE["inverse_trig"] = BenchmarkGroup()
const m1_small = Multicomplex(0.5, 0.3)
const m2_small = Multicomplex(0.5, 0.3, 0.2, 0.1)
SUITE["inverse_trig"]["asin N=1"] = @benchmarkable asin($m1_small)
SUITE["inverse_trig"]["asin N=2"] = @benchmarkable asin($m2_small)
SUITE["inverse_trig"]["acos N=1"] = @benchmarkable acos($m1_small)
SUITE["inverse_trig"]["acos N=2"] = @benchmarkable acos($m2_small)
SUITE["inverse_trig"]["atan N=1"] = @benchmarkable atan($m1_small)
SUITE["inverse_trig"]["atan N=2"] = @benchmarkable atan($m2_small)

# Inverse hyperbolic functions
SUITE["inverse_hyperbolic"] = BenchmarkGroup()
SUITE["inverse_hyperbolic"]["asinh N=1"] = @benchmarkable asinh($m1_small)
SUITE["inverse_hyperbolic"]["asinh N=2"] = @benchmarkable asinh($m2_small)
SUITE["inverse_hyperbolic"]["atanh N=1"] = @benchmarkable atanh($m1_small)
SUITE["inverse_hyperbolic"]["atanh N=2"] = @benchmarkable atanh($m2_small)

# Fold and isabient
SUITE["fold"] = BenchmarkGroup()
SUITE["fold"]["N=1"] = @benchmarkable fold($m1)
SUITE["fold"]["N=2"] = @benchmarkable fold($m2)
SUITE["fold"]["N=3"] = @benchmarkable fold($m3)

const m_abient = 1.0 + im1*im2
SUITE["fold"]["isabient true"] = @benchmarkable isabient($m_abient)
SUITE["fold"]["isabient false"] = @benchmarkable isabient($m2)

# Conjugation
SUITE["conjugation"] = BenchmarkGroup()
SUITE["conjugation"]["N=1"] = @benchmarkable conj($m1)
SUITE["conjugation"]["N=2"] = @benchmarkable conj($m2)
SUITE["conjugation"]["N=3"] = @benchmarkable conj($m3)

# Absolute value
SUITE["abs"] = BenchmarkGroup()
SUITE["abs"]["abs N=1"] = @benchmarkable abs($m1)
SUITE["abs"]["abs N=2"] = @benchmarkable abs($m2)
SUITE["abs"]["abs N=3"] = @benchmarkable abs($m3)
SUITE["abs"]["abs2 N=1"] = @benchmarkable abs2($m1)
SUITE["abs"]["abs2 N=2"] = @benchmarkable abs2($m2)
SUITE["abs"]["abs2 N=3"] = @benchmarkable abs2($m3)

# Component extraction
SUITE["components"] = BenchmarkGroup()
SUITE["components"]["real N=1"] = @benchmarkable real($m1)
SUITE["components"]["real N=2"] = @benchmarkable real($m2)
SUITE["components"]["real N=3"] = @benchmarkable real($m3)
SUITE["components"]["imag N=1"] = @benchmarkable imag($m1)
SUITE["components"]["imag N=2"] = @benchmarkable imag($m2)
SUITE["components"]["imag N=3"] = @benchmarkable imag($m3)
SUITE["components"]["realest N=2"] = @benchmarkable realest($m2)
SUITE["components"]["realest N=3"] = @benchmarkable realest($m3)
SUITE["components"]["flat N=2"] = @benchmarkable flat($m2)
SUITE["components"]["flat N=3"] = @benchmarkable flat($m3)

# Matrix representation (expensive but fundamental)
SUITE["matrep"] = BenchmarkGroup()
SUITE["matrep"]["N=0"] = @benchmarkable matrep($m0)
SUITE["matrep"]["N=1"] = @benchmarkable matrep($m1)
SUITE["matrep"]["N=2"] = @benchmarkable matrep($m2)
SUITE["matrep"]["N=3"] = @benchmarkable matrep($m3)
SUITE["matrep"]["N=4"] = @benchmarkable matrep($m4)

# Random number generation
SUITE["random"] = BenchmarkGroup()
SUITE["random"]["rand N=1"] = @benchmarkable rand(Multicomplex{Float64,1,2})
SUITE["random"]["rand N=2"] = @benchmarkable rand(Multicomplex{Float64,2,4})
SUITE["random"]["rand N=3"] = @benchmarkable rand(Multicomplex{Float64,3,8})
SUITE["random"]["randn N=1"] = @benchmarkable randn(Multicomplex{Float64,1,2})
SUITE["random"]["randn N=2"] = @benchmarkable randn(Multicomplex{Float64,2,4})
SUITE["random"]["randn N=3"] = @benchmarkable randn(Multicomplex{Float64,3,8})

# Array operations
SUITE["arrays"] = BenchmarkGroup()
const arr1 = [Multicomplex(1.0, 2.0) for _ in 1:100]
const arr2 = [Multicomplex(1.0, 2.0, 3.0, 4.0) for _ in 1:100]
SUITE["arrays"]["sum N=1 array"] = @benchmarkable sum($arr1)
SUITE["arrays"]["sum N=2 array"] = @benchmarkable sum($arr2)
SUITE["arrays"]["broadcast add N=1"] = @benchmarkable $arr1 .+ $arr1
SUITE["arrays"]["broadcast add N=2"] = @benchmarkable $arr2 .+ $arr2
SUITE["arrays"]["broadcast mul scalar N=1"] = @benchmarkable 2.0 .* $arr1
SUITE["arrays"]["broadcast mul scalar N=2"] = @benchmarkable 2.0 .* $arr2
