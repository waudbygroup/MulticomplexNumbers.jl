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
