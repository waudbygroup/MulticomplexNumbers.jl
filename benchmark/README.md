# Benchmarking MulticomplexNumbers.jl

This directory contains performance benchmarks for the MulticomplexNumbers.jl package using [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl).

## Setup

First, install BenchmarkTools:

```julia
using Pkg
Pkg.add("BenchmarkTools")
```

## Running Benchmarks

### Run all benchmarks

```julia
using BenchmarkTools
using MulticomplexNumbers

include("benchmark/benchmarks.jl")
results = run(SUITE)
```

### Run specific benchmark groups

```julia
# Run only arithmetic benchmarks
results = run(SUITE["arithmetic"])

# Run only multiplication benchmarks
results = run(SUITE["multiplication"])

# Run only transcendental function benchmarks
results = run(SUITE["transcendental"])
```

### View results

```julia
# Show results in a readable format
using BenchmarkTools
show(stdout, MIME("text/plain"), results)

# Get median times
median(results)

# Get minimum times
minimum(results)

# Compare two sets of results
judge(results_new, results_old)
```

## Benchmark Categories

The benchmark suite is organized into the following categories:

- **constructors**: Different ways to create multicomplex numbers
- **arithmetic**: Basic operations (+, -, *, /, scalar operations)
- **multiplication**: Multicomplex-multicomplex multiplication for different orders
- **division**: Multicomplex division for different orders
- **powers**: Integer and real powers
- **transcendental**: exp, log, sqrt
- **trigonometric**: sin, cos, tan, sec, csc, cot
- **hyperbolic**: sinh, cosh, tanh, sech, csch, coth
- **inverse_trig**: asin, acos, atan, asec, acsc, acot
- **inverse_hyperbolic**: asinh, acosh, atanh, asech, acsch, acoth
- **fold**: fold operator and isabient function
- **conjugation**: Complex conjugation
- **abs**: Absolute value and abs2
- **components**: Component extraction (real, imag, realest, flat)
- **matrep**: Matrix representation generation
- **random**: Random number generation (rand, randn)
- **arrays**: Array operations and broadcasting

## Performance Tips

Based on the benchmarks, here are some performance considerations:

1. **Matrix operations are expensive**: Operations like multiplication, division, powers, and transcendental functions use the matrix representation, which is computationally expensive for higher orders.

2. **Order scaling**: Performance scales with the order N:
   - N=0: Very fast (essentially scalar operations)
   - N=1: Fast (2 components, similar to Complex)
   - N=2: Moderate (4 components, 4×4 matrix operations)
   - N=3: Slower (8 components, 8×8 matrix operations)
   - N=4+: Much slower (exponential growth in matrix size)

3. **Cheap operations**: Component-wise operations (addition, subtraction, scalar multiplication) are very fast and scale linearly with the number of components.

4. **Expensive operations**: Matrix-based operations (multiplication, division, transcendental functions) scale poorly with order due to matrix operations.

## Saving and Loading Results

To track performance over time:

```julia
using BenchmarkTools, JSON

# Run benchmarks
results = run(SUITE)

# Save to file
BenchmarkTools.save("benchmark_results.json", results)

# Load from file
results_old = BenchmarkTools.load("benchmark_results.json")[1]

# Compare with new results
results_new = run(SUITE)
comparison = judge(results_new, results_old)
```

## Continuous Benchmarking

For CI/CD integration, consider using [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl):

```julia
using PkgBenchmark

# Benchmark current state
results = benchmarkpkg("MulticomplexNumbers")

# Compare against a specific commit or tag
judge("MulticomplexNumbers", "v0.1.0")
```

## Contributing

When adding new features to MulticomplexNumbers.jl:

1. Add corresponding benchmarks to `benchmarks.jl`
2. Run benchmarks before and after changes
3. Document any significant performance impacts
4. Consider adding performance regression tests for critical operations
