# CLAUDE.md - AI Assistant Guide for MulticomplexNumbers.jl

This document provides comprehensive guidance for AI assistants working with the MulticomplexNumbers.jl codebase.

## Table of Contents

1. [Repository Overview](#repository-overview)
2. [Project Structure](#project-structure)
3. [Core Concepts](#core-concepts)
4. [Key Files Reference](#key-files-reference)
5. [Development Workflow](#development-workflow)
6. [Testing Strategy](#testing-strategy)
7. [Code Conventions](#code-conventions)
8. [CI/CD Pipeline](#cicd-pipeline)
9. [Documentation](#documentation)
10. [Common Development Tasks](#common-development-tasks)
11. [Important Notes for AI Assistants](#important-notes-for-ai-assistants)

---

## Repository Overview

**Package Name**: MulticomplexNumbers.jl
**Purpose**: A Julia package for representing multicomplex numbers and performing multicomplex algebra
**Author**: Chris Waudby (c.waudby@ucl.ac.uk)
**License**: MIT
**Julia Version**: 1.9+

### What are Multicomplex Numbers?

Multicomplex numbers are a generalization of complex numbers, recursively defined to contain multiple imaginary numbers (i₁, i₂, etc.). Unlike Clifford algebras, these numbers commute: i₁i₂ = i₂i₁.

### Key Dependencies

- **LinearAlgebra**: Standard Julia linear algebra
- **StaticArrays**: For efficient fixed-size arrays (core to performance)
- **PrecompileTools**: For precompilation optimization
- **FFTW** (weak dependency): FFT support via package extension

---

## Project Structure

```
MulticomplexNumbers.jl/
├── src/                          # Source code
│   ├── MulticomplexNumbers.jl    # Main module file (exports, includes)
│   ├── base.jl                   # Core type definition, constructors, conversions
│   ├── arithmetic.jl             # Arithmetic operations (+, -, *, /, ^, exp, log, sqrt)
│   ├── representations.jl        # Matrix representations and ascomplex conversions
│   ├── io.jl                     # Display and printing
│   └── precompile.jl            # Precompilation statements
├── ext/                          # Package extensions (Julia 1.9+)
│   └── FFTWExt.jl               # FFTW integration for FFT operations
├── test/                         # Test suite
│   ├── runtests.jl              # Main test runner (uses SafeTestsets)
│   ├── base_test.jl             # Tests for core functionality
│   └── fftwext_test.jl          # Tests for FFTW extension
├── docs/                         # Documentation
│   ├── make.jl                  # Documenter.jl build script
│   ├── Project.toml             # Documentation dependencies
│   └── src/                     # Documentation source files
│       ├── index.md
│       ├── background.md
│       ├── userguide.md
│       ├── examples.md
│       └── api.md
├── .github/workflows/           # CI/CD
│   ├── Runtests.yml            # Test automation
│   ├── Documenter.yml          # Documentation deployment
│   ├── CompatHelper.yml        # Dependency updates
│   └── TagBot.yml              # Release automation
├── Project.toml                 # Package metadata and dependencies
├── README.md                    # User-facing readme
└── LICENSE                      # MIT license
```

---

## Core Concepts

### The Multicomplex Type

```julia
struct Multicomplex{T,N,C} <: Number
    value::SVector{C,T}
end
```

**Type Parameters**:
- `T`: Base scalar type (must be `<:Real`, checked by `can_multicomplex()`)
- `N`: Order of the multicomplex number (number of imaginary units)
- `C`: Number of components (always C = 2^N)

**Storage**: Uses `SVector` from StaticArrays for efficient, stack-allocated storage.

### Imaginary Units

Pre-defined constants for convenience:
- `im1`: Multicomplex{1} - equivalent to standard `im`
- `im2`: Multicomplex{2} - second imaginary unit
- `im3` through `im6`: Higher order units

### Component Ordering

Components are stored in a specific order. For example, Multicomplex{2} has 4 components:
- `[1]`: real part (no imaginary units)
- `[2]`: i₁ coefficient
- `[3]`: i₂ coefficient
- `[4]`: i₁i₂ coefficient

### Matrix Representation

Multicomplex numbers can be represented as matrices, which is used for:
- Multiplication
- Division
- Powers
- Exponentials, logarithms, square roots

The `matrep()` function converts a multicomplex number to its matrix form.

---

## Key Files Reference

### src/MulticomplexNumbers.jl (24 lines)

**Purpose**: Main module definition and exports

**Exports**:
- `Multicomplex`: The main type
- `im1` through `im6`: Imaginary units
- `order`: Get the order N of a multicomplex number
- `flat`: Get components as a vector
- `component`: Extract a specific component
- `realest`: Get the most real component
- `matrep`: Get matrix representation
- `ascomplex`: View as complex array

**Structure**: Includes other source files in order:
1. `base.jl` - type definitions
2. `io.jl` - display methods
3. `representations.jl` - matrix and complex conversions
4. `arithmetic.jl` - mathematical operations
5. `precompile.jl` - precompilation directives

### src/base.jl (195 lines)

**Purpose**: Core type definition, constructors, conversions, and basic operations

**Key Functions**:
- `can_multicomplex(::Type)`: Type trait for valid scalar types
- Multiple constructor methods for different input types
- Conversion and promotion rules
- Component extraction (`real()`, `imag()`, `realest()`, `component()`)
- Testing functions (`isreal()`, `isinteger()`, `isfinite()`, etc.)
- `norm()`: Euclidean norm of components

**Important Pattern**: Uses `@boundscheck` for performance-critical code

### src/arithmetic.jl (144 lines)

**Purpose**: All arithmetic operations

**Operations Implemented**:
- Comparison: `==`, `isequal`, `hash`
- Conjugation: `conj()`
- Absolute values: `abs()`, `abs2()`
- Basic arithmetic: `+`, `-`, `*`, `/` (scalar and multicomplex)
- Powers: `^` (via matrix representation)
- Exponentials: `exp()`
- Logarithms: `log()`
- Square root: `sqrt()`

**Performance Note**: Most operations use the matrix representation, which is computed on-the-fly.

### src/representations.jl (334 lines)

**Purpose**: Matrix representations and complex array conversions

**Key Functions**:
- `matrep(m)`: Convert multicomplex to matrix representation
  - Explicitly defined for N=0,1,2,3
  - Recursive definition for higher orders
- `ascomplex(A, unit)`: View multicomplex array as complex array
  - Maps i(unit) → im
  - Returns a view (no copying)
  - Handles dimension permutations
- `unsafe_ascomplex!()` and `unsafe_fromcomplex!()`: In-place conversions for FFTW

### src/io.jl (58 lines)

**Purpose**: Display and printing methods

**Functionality**:
- Custom `show()` methods for REPL display
- Pretty-printing with imaginary unit subscripts (i₁, i₂, etc.)
- Handles different output contexts (inline vs multiline)

### ext/FFTWExt.jl (48 lines)

**Purpose**: FFTW integration via package extension (Julia 1.9+)

**Key Function**: `fft!(A, unit, dims)`
- Performs in-place FFT on multicomplex arrays
- Uses `unsafe_ascomplex!()` to convert to complex, perform FFT, then convert back
- Handles proper dimension mapping for different multicomplex orders

**Design Pattern**: Package extensions only load when FFTW is explicitly loaded by the user.

---

## Development Workflow

### Setting Up Development Environment

```julia
# Clone the repository
git clone https://github.com/waudbygroup/MulticomplexNumbers.jl.git
cd MulticomplexNumbers.jl

# Activate the project
julia --project=.

# In Julia REPL:
]  # Enter package mode
instantiate  # Install dependencies
```

### Working with the Package

```julia
# Activate and develop the package
]activate .
]instantiate

# Run tests
]test

# Build documentation
]activate docs
include("docs/make.jl")
```

### Adding Features

1. **Determine the right file**:
   - New type functionality → `src/base.jl`
   - Mathematical operations → `src/arithmetic.jl`
   - Matrix/representation work → `src/representations.jl`
   - Display changes → `src/io.jl`

2. **Add tests**: Every new feature must have corresponding tests in `test/base_test.jl`

3. **Update exports**: Add new public functions to exports in `src/MulticomplexNumbers.jl`

4. **Document**: Add docstrings to public functions

---

## Testing Strategy

### Test Structure

```julia
# test/runtests.jl - main entry point
using MulticomplexNumbers
using Test
using SafeTestsets

@safetestset "MulticomplexNumbers" begin
    include("base_test.jl")
end

@safetestset "FFTWExt" begin
    include("fftwext_test.jl")
end
```

**SafeTestsets**: Each test file runs in its own module to prevent namespace pollution.

### Running Tests

```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Or in Julia REPL
]test
```

### Test Coverage

- Tests run on macOS-latest (configurable in CI)
- Coverage reports uploaded to Codecov
- Currently testing Julia 1.9 (nightly commented out)

### Testing FFTW Extension

```julia
# The extension is automatically loaded when FFTW is loaded
]activate --temp
]dev .
]add FFTW
using MulticomplexNumbers, FFTW
# Now fft! methods are available
```

---

## Code Conventions

### Style Guidelines

1. **Type Stability**: All functions should be type-stable for performance
2. **StaticArrays**: Use `SVector` and `SMatrix` for fixed-size arrays
3. **Docstrings**: Use triple-quoted docstrings for public functions
4. **Mathematical Notation**: Use Julia's unicode support (e.g., ℂ, subscripts)

### Naming Conventions

- **Types**: PascalCase (e.g., `Multicomplex`)
- **Functions**: snake_case for internal, or descriptive names (e.g., `matrep`, `can_multicomplex`)
- **Constants**: lowercase for mathematical constants (e.g., `im1`, `im2`)
- **Type parameters**: Single uppercase letters (T, N, C, S)

### Performance Patterns

```julia
# Use @boundscheck for performance-critical code
@boundscheck begin
    # bounds checking code
end

# Use type parameters to specialize functions
function foo(m::Multicomplex{T,N,C}) where {T,N,C}
    # compiler can optimize based on N and C
end

# Prefer StaticArrays for fixed-size operations
SVector{4,T}(a, b, c, d)  # Better than [a,b,c,d]
```

### Error Handling

```julia
# Use specific exceptions
throw(ArgumentError("Cannot create..."))
throw(InexactError(nameof(T), T, m))

# Use @noinline for error paths (performance)
@noinline function throw_cannot_multicomplex(T::Type)
    throw(ArgumentError("..."))
end
```

---

## CI/CD Pipeline

### GitHub Actions Workflows

#### Runtests.yml
**Trigger**: Push and pull requests
**Purpose**: Run test suite
**Matrix**:
- Julia versions: 1.9 (nightly commented out)
- OS: macOS-latest (ubuntu-latest and windows-latest commented out)
- Architecture: x64

**Steps**:
1. Checkout code
2. Setup Julia
3. Build package
4. Run tests with annotations
5. Process coverage
6. Upload to Codecov

#### Documenter.yml
**Trigger**: Push to main, tags, pull requests
**Purpose**: Build and deploy documentation to GitHub Pages

#### CompatHelper.yml
**Purpose**: Automated dependency compatibility updates
**Frequency**: Scheduled (creates PRs for dependency updates)

#### TagBot.yml
**Purpose**: Automatic GitHub releases when package is tagged
**Integration**: Works with Julia Registry

### Code Coverage

- Coverage tracking via Codecov
- Token: V9ND8Y3R8A (stored in repository secrets)
- Reports available at: https://codecov.io/gh/waudbygroup/MulticomplexNumbers.jl

---

## Documentation

### Documentation Structure

Built with Documenter.jl, hosted on GitHub Pages:
- **Stable**: https://waudbygroup.github.io/MulticomplexNumbers.jl/stable
- **Dev**: https://waudbygroup.github.io/MulticomplexNumbers.jl/dev

### Documentation Files

- `index.md`: Landing page, installation, citations
- `background.md`: Mathematical background on multicomplex numbers
- `userguide.md`: Usage guide and tutorials
- `examples.md`: Example code and applications
- `api.md`: API reference (auto-generated from docstrings)

### Building Documentation Locally

```bash
cd docs
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. make.jl
```

Output will be in `docs/build/`.

### Adding Documentation

1. Add docstrings to public functions:
   ```julia
   """
       function_name(args)

   Brief description.

   # Arguments
   - `arg1`: Description

   # Returns
   Description of return value

   # Examples
   ```julia
   julia> example_code
   ```
   """
   function function_name(args)
   ```

2. Update relevant .md files in `docs/src/`

3. Rebuild documentation to verify

---

## Common Development Tasks

### Adding a New Arithmetic Operation

```julia
# 1. Add to src/arithmetic.jl
"""
    my_operation(m::Multicomplex)

Description of the operation.
"""
function Base.my_operation(m::Multicomplex{T,N,C}) where {T,N,C}
    # Implementation
    Multicomplex{N}(result_vector)
end

# 2. Add tests to test/base_test.jl
@testset "my_operation" begin
    m = Multicomplex(1.0, 2.0)
    result = my_operation(m)
    @test result ≈ expected_value
end

# 3. Export if public (src/MulticomplexNumbers.jl)
export my_operation

# 4. Add documentation (docs/src/api.md)
```

### Adding a New Constructor

```julia
# Add to src/base.jl in the constructor section
Multicomplex(special_input_type) = Multicomplex{N}(SVector(...))

# Add promotion rules if needed
Base.promote_rule(::Type{Multicomplex{T,N,C}}, ::Type{NewType}) where {T,N,C} = ...
```

### Working with the FFTW Extension

```julia
# Extension code is in ext/FFTWExt.jl
# It only loads when FFTW is loaded

# To test:
using MulticomplexNumbers
using FFTW  # This triggers the extension load

# Now fft! is available for Multicomplex arrays
A = [Multicomplex(1.0, 2.0) for _ in 1:10]
fft!(A, 1)  # FFT along unit 1
```

### Debugging Performance Issues

```julia
# Use @code_warntype to check type stability
@code_warntype my_function(args)

# Use BenchmarkTools for timing
using BenchmarkTools
@benchmark my_function($args)

# Profile code
using Profile
@profile my_function(args)
Profile.print()
```

---

## Important Notes for AI Assistants

### Type System Constraints

1. **Scalar Type Restrictions**: Only `<:Real` types can be used as the base type `T`. This is enforced by `can_multicomplex(::Type)`.

2. **Component Count**: Always C = 2^N. This is checked in the constructor.

3. **StaticArrays Required**: The `value` field must be an `SVector{C,T}`. Never use regular arrays.

### Performance Critical Paths

1. **Matrix Representation**: `matrep()` is called frequently in arithmetic operations. It's expensive for high orders.

2. **Avoid Copying**: Use views and in-place operations where possible. The `ascomplex()` function returns a view, not a copy.

3. **Type Stability**: All operations should infer concrete types. Use `@code_warntype` to verify.

### Common Pitfalls

1. **Component Indexing**: Components are 1-indexed and follow a specific order. Don't assume arbitrary ordering.

2. **Order Mismatches**: Operations between different orders require promotion. The type system handles this automatically.

3. **Complex vs Multicomplex**: `Multicomplex{1}` is similar to `Complex` but not identical. They can be converted but are different types.

4. **Extension Loading**: FFTW methods only exist when FFTW is loaded. Check for extension availability before using FFTW-specific features.

### Testing Requirements

1. **SafeTestsets**: Always use `@safetestset` for new test groups to prevent namespace pollution.

2. **Numerical Accuracy**: Use `≈` (isapprox) for floating-point comparisons, not `==`.

3. **Edge Cases**: Test with:
   - Zero values
   - Integer types
   - High precision types (BigFloat)
   - Different orders (N=0,1,2,3,4)

### Git Workflow

1. **Branching**: Development happens on feature branches starting with `claude/`

2. **Commits**: Should be descriptive and focused. One logical change per commit.

3. **Push**: Always use `git push -u origin <branch-name>` for new branches.

4. **CI Must Pass**: All tests must pass before merging. Check the GitHub Actions status.

### Documentation Standards

1. **Public API**: All exported functions must have docstrings.

2. **Examples**: Include at least one example in docstrings for non-trivial functions.

3. **Math Notation**: Use LaTeX in docstrings for mathematical formulas (wrapped in `` `` or ``` ``` for display math).

### When Adding Dependencies

1. **Weak Dependencies**: Use weak dependencies (`[weakdeps]`) for optional features.

2. **Extensions**: Use package extensions (Julia 1.9+) instead of Requires.jl.

3. **Version Constraints**: Add appropriate `[compat]` entries in Project.toml.

4. **Documentation Dependencies**: Keep docs/Project.toml separate from main Project.toml.

### Known TODOs and Future Work

From README.md:
- ~~FFT support~~ (DONE - implemented via FFTWExt)

Potential future enhancements:
- Support for more Julia versions in CI
- Support for more operating systems in CI (currently macOS only)
- Additional mathematical functions
- Performance optimizations for high-order multicomplex numbers

### References and Resources

**External Documentation**:
- NIST report on multicomplex algebra: https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf
- NIST C++ implementation: https://github.com/usnistgov/multicomplex
- Academic paper (Casado & Hewson, 2020): http://dx.doi.org/10.1145/3378542

**Julia Development Resources**:
- Package creation: https://jaantollander.com/post/how-to-create-software-packages-with-julia-language/
- Documenter.jl: https://juliadocs.github.io/Documenter.jl/stable/
- Package tutorial: https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/11-developing-julia-packages

### Quick Reference Commands

```bash
# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Build docs locally
julia --project=docs docs/make.jl

# Start REPL with package
julia --project=.

# Check code coverage locally
julia --project=. --code-coverage=user -e 'using Pkg; Pkg.test()'

# Format code (if JuliaFormatter is available)
using JuliaFormatter; format("src/")
```

---

## Summary

MulticomplexNumbers.jl is a well-structured, performance-oriented Julia package implementing multicomplex number algebra. The codebase prioritizes:

1. **Performance**: Heavy use of StaticArrays, type stability, and careful memory management
2. **Correctness**: Comprehensive test suite with SafeTestsets
3. **Maintainability**: Clear code organization, good documentation, automated CI/CD
4. **Extensibility**: Package extension system for optional dependencies

When working with this codebase, always consider type stability, maintain the test suite, and follow the established patterns for type parameters and StaticArrays usage.

---

## Codebase Analysis & Critique

### Overall Assessment

**Grade: A- (87/100)**

This is a **high-quality, well-designed package** with excellent fundamentals. The type system is clever, the code is clean, and the architecture is sound.

**Code Quality Metrics**:
- **Lines of Code**: ~1,221 (compact, well-organized)
- **Test Coverage**: High (good use of `@inferred`)
- **Type Stability**: ✅ Excellent (all operations type-stable)
- **Documentation**: ⚠️ Moderate (API docs good, tutorials needed)
- **Performance**: ✅ Good (StaticArrays, minimal allocations)
- **Modularity**: ✅ Excellent (clean file organization)
- **CI/CD**: ⚠️ Limited (only macOS tested)
- **Dependencies**: ✅ Minimal (good use of weak deps)

### Strengths

1. **Excellent Type Design**
   - Clever use of StaticArrays for performance (stack allocation)
   - Type parameters (T, N, C) well-chosen and enforce invariants at compile time
   - Proper Julia Number interface implementation
   - Good use of promotion rules and conversion methods

2. **Strong Architecture**
   - Clean separation of concerns across files (base, arithmetic, representations, io)
   - Matrix representation approach is mathematically sound
   - Package extension for FFTW is modern Julia 1.9+ best practice
   - Recursive `matrep()` definition is elegant

3. **Good Test Coverage**
   - Tests use `@inferred` to verify type stability
   - SafeTestsets prevent namespace pollution
   - Good coverage of edge cases (N=0,1,2,3,4)
   - FFTW extension has comprehensive tests

4. **Performance Awareness**
   - `@boundscheck` used appropriately
   - Views instead of copies (`ascomplex` returns views)
   - Precompilation workload defined
   - Type-stable operations throughout

### Weaknesses & Issues

#### 1. Limited CI Coverage
Currently only testing on macOS-latest with Julia 1.9. No testing on Linux/Windows, only single Julia version.

**Location**: `.github/workflows/Runtests.yml:10`

#### 2. Missing Mathematical Functions
The package implements only basic operations:
- ✅ Has: `+, -, *, /, ^, exp, log, sqrt, abs, conj`
- ❌ Missing: `sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, atan2`
- ❌ Missing: Special functions that other Number types support

**Impact**: Users cannot use multicomplex numbers in contexts requiring these functions.

#### 3. Incomplete FFTW Extension
Only supports N=1,2,3. Higher orders (N≥4) not supported for FFT.

**Location**: `ext/FFTWExt.jl:39`

#### 4. No Performance Benchmarks
- No benchmark suite to track performance regressions
- No comparison with other implementations (NIST C++, Matlab)
- `matrep()` is expensive but not profiled/optimized

#### 5. Documentation Gaps
- Missing examples of practical applications
- No performance guidelines for users
- Mathematical background could be more detailed
- No migration guide from Complex to Multicomplex

#### 6. Limited Error Messages
Error messages could be more helpful with suggestions.

**Example**: `throw(ArgumentError("unsupported multicomplex order"))` could suggest filing an issue.

#### 7. No Broadcasting Optimizations
Standard broadcasting works but there may be optimization opportunities for multicomplex-specific operations.

#### 8. Missing Utility Functions
- No `zero(::Type{Multicomplex{T,N,C}})` or `one(...)` constructors
- No `rand(Multicomplex{N})` for random numbers
- No conversion to/from other formats (JSON, HDF5, etc.)

#### 9. Test Organization
All tests in one 346-line file (`base_test.jl`) makes it harder to navigate.

#### 10. README TODO Still Present
The README.md still contains "TODO: FFT" even though FFT support has been implemented via FFTWExt!

**Location**: `README.md:22`

---

## Recommended Next Steps

### Priority 1: High Impact, Easy Wins

#### 1.1 Update README.md
- ✅ Remove "TODO: FFT" (it's done!)
- Add installation instructions
- Add quick start example
- Add badges for docs/CI/coverage
- Update references to show FFT is implemented

**Estimated effort**: 30 minutes

#### 1.2 Expand CI Matrix
Add testing on multiple platforms and Julia versions:
```yaml
os: [ubuntu-latest, windows-latest, macOS-latest]
julia-version: ['1.9', '1.10', 'nightly']
```

**Location**: `.github/workflows/Runtests.yml`
**Estimated effort**: 15 minutes

#### 1.3 Add Missing Base Functions
```julia
# In src/base.jl
Base.zero(::Type{Multicomplex{T,N,C}}) where {T,N,C} =
    Multicomplex{N}(zeros(SVector{C,T}))
Base.one(::Type{Multicomplex{T,N,C}}) where {T,N,C} =
    Multicomplex{N}(SVector(one(T), zeros(T, C-1)...))
```

**Estimated effort**: 20 minutes

#### 1.4 Split Test Files
Reorganize tests into separate files:
```
test/
├── runtests.jl
├── constructors_test.jl
├── arithmetic_test.jl
├── representations_test.jl
├── io_test.jl
└── fftwext_test.jl
```

**Estimated effort**: 30 minutes

#### 1.5 Improve Error Messages
Make error messages more helpful:
```julia
throw(ArgumentError(
    "FFT for multicomplex order $N is not yet supported. " *
    "Currently supported orders: 1, 2, 3. " *
    "Please file an issue at https://github.com/waudbygroup/MulticomplexNumbers.jl/issues"
))
```

**Estimated effort**: 20 minutes

#### 1.6 Extend FFTW Support to N=4,5
Add support for higher order multicomplex numbers in FFT operations:
```julia
# In ext/FFTWExt.jl
elseif N == 4
    if unit == 1
        fft!(w, dims)
    elseif unit == 2
        d = dims .+ 1
        fft!(w, d)
    elseif unit == 3
        d = dims .+ 2
        fft!(w, d)
    else
        d = dims .+ 3
        fft!(w, d)
    end
# Similar pattern for N == 5
```

**Location**: `ext/FFTWExt.jl:39`
**Estimated effort**: 3 hours

**Rationale**: FFT is a core feature mentioned in the package description. Completing support for higher orders removes a significant limitation and improves the package's utility for scientific computing applications.

### Priority 2: Mathematical Completeness

#### 2.1 Implement Trigonometric Functions
Add via matrix representation (same pattern as exp, log, sqrt):
```julia
Base.sin(m::Multicomplex{T,N,C}) where {T,N,C} =
    Multicomplex{N}(sin(matrep(m))[SVector{C}(SOneTo(C))])
Base.cos(m::Multicomplex{T,N,C}) where {T,N,C} =
    Multicomplex{N}(cos(matrep(m))[SVector{C}(SOneTo(C))])
# Similarly for tan, cot, sec, csc
```

**Location**: Create new file `src/trigonometric.jl`
**Estimated effort**: 2 hours

#### 2.2 Add Hyperbolic Functions
```julia
Base.sinh, Base.cosh, Base.tanh, Base.asinh, Base.acosh, Base.atanh
```

**Estimated effort**: 1 hour

#### 2.3 Implement Inverse Trig Functions
```julia
Base.asin, Base.acos, Base.atan
```

**Estimated effort**: 1 hour

### Priority 3: Performance & Quality

#### 3.1 Add Benchmark Suite
Create `benchmark/benchmarks.jl`:
```julia
using BenchmarkTools
using MulticomplexNumbers

SUITE = BenchmarkGroup()
SUITE["arithmetic"] = BenchmarkGroup()
SUITE["arithmetic"]["multiplication"] = @benchmarkable m1 * m2 setup=(m1=Multicomplex(1.0,2.0); m2=Multicomplex(3.0,4.0))
# Compare against Complex for N=1
```

**Estimated effort**: 3 hours

#### 3.2 Profile and Optimize `matrep()`
- Consider caching for repeated operations
- Benchmark explicit vs recursive definitions
- Document performance characteristics
- Add benchmarks to track regressions

**Estimated effort**: 4 hours

#### 3.3 Add Type Piracy Checks
```julia
# In test/aqua_test.jl
using Aqua
@testset "Aqua.jl" begin
    Aqua.test_all(MulticomplexNumbers)
end
```

**Dependencies**: Add Aqua.jl to test dependencies
**Estimated effort**: 1 hour

### Priority 4: Enhanced Usability

#### 4.1 Add Random Number Generation
```julia
# In src/base.jl
using Random
Base.rand(rng::AbstractRNG, ::Random.SamplerType{Multicomplex{T,N,C}}) where {T,N,C} =
    Multicomplex{N}(rand(rng, SVector{C,T}))
```

**Estimated effort**: 30 minutes

#### 4.2 Add Conversion Utilities
```julia
# To/from nested Complex numbers
to_nested_complex(m::Multicomplex{T,2}) = Complex(real(m), imag(m))
from_nested_complex(z::Complex{<:Complex}) = Multicomplex(real(real(z)), imag(real(z)), real(imag(z)), imag(imag(z)))
```

**Estimated effort**: 1 hour

#### 4.3 Broadcasting Optimizations
Specialized broadcasting for common patterns to avoid unnecessary allocations.

**Estimated effort**: 3 hours

### Priority 5: Documentation

#### 5.1 Expand Documentation
- Add "Performance Tips" page
- Add "Migration from Complex" guide
- Add "Applications" page with real examples (numerical differentiation!)
- Add inline examples to all docstrings

**Estimated effort**: 4 hours

#### 5.2 Add Tutorials
Create tutorial notebooks:
- Numerical differentiation (the main use case!)
- Signal processing with multicomplex FFT
- Solving PDEs with multicomplex numbers

**Estimated effort**: 6 hours

#### 5.3 Create Comparison Table
Add to documentation:
```markdown
| Feature | Complex | Multicomplex{1} | Multicomplex{2} |
|---------|---------|-----------------|-----------------|
| Storage | 16 bytes| 16 bytes        | 32 bytes        |
| # Imaginary units | 1 | 1 | 2 |
| FFT support | Yes | Yes | Yes |
```

**Estimated effort**: 1 hour

### Priority 6: Ecosystem Integration

#### 6.1 Add ChainRules Integration
For automatic differentiation compatibility:
```julia
# Create ext/ChainRulesExt.jl
using ChainRulesCore
# Define frules and rrules for key operations
```

**Dependencies**: Add ChainRules to weakdeps
**Estimated effort**: 4 hours

#### 6.2 Add StructTypes Integration
For JSON serialization:
```julia
# Create ext/StructTypesExt.jl
using StructTypes
StructTypes.StructType(::Type{<:Multicomplex}) = StructTypes.Struct()
```

**Estimated effort**: 1 hour

---

## Strategic Roadmap

### Short Term (1-2 weeks)
1. Fix README and update TODO list
2. Expand CI to Linux/Windows
3. Add `zero`, `one`, and basic utilities
4. Split test files for maintainability
5. Improve error messages
6. Extend FFTW support to N=4,5

**Total effort**: ~7 hours

### Medium Term (1-2 months)
1. Implement trigonometric functions
2. Implement hyperbolic and inverse trig functions
3. Add comprehensive benchmarks
4. Write performance optimization guide
5. Add ChainRules for AD compatibility
6. Create tutorial notebooks

**Total effort**: ~29 hours

### Long Term (3-6 months)
1. Performance optimization campaign
2. Write academic paper on implementation
3. Create comparison with NIST/Matlab implementations
4. Add GPU support (via CUDA.jl extension)
5. Integrate with DifferentialEquations.jl ecosystem

**Total effort**: ~60+ hours

---

## Development Checklist for Contributors

When adding new features, ensure you:

- [ ] Add comprehensive tests with `@inferred` checks
- [ ] Add docstrings with examples
- [ ] Update exports in `src/MulticomplexNumbers.jl`
- [ ] Update CLAUDE.md if adding significant functionality
- [ ] Run tests on Linux, macOS, Windows (or wait for CI)
- [ ] Check type stability with `@code_warntype`
- [ ] Add to precompile workload if commonly used
- [ ] Update documentation in `docs/src/`
- [ ] Add benchmarks if performance-critical
- [ ] Follow existing code style and patterns

---

## Package Maturity Assessment

**Recommendation**: This package is ready for **v1.0 release** after addressing Priority 1 items.

The codebase demonstrates:
- ✅ Strong Julia programming skills
- ✅ Mathematical sophistication
- ✅ Good software engineering practices
- ✅ Active maintenance (recent commits)
- ⚠️ Limited real-world usage (needs promotion)
- ⚠️ Some missing standard features

**Next milestone**: Address Priority 1 and Priority 2 items, then release v1.0.0.
