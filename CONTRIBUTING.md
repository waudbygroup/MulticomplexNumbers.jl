# Contributing to MulticomplexNumbers.jl

Thanks for your interest in improving MulticomplexNumbers.jl! Contributions of
all kinds are welcome — bug reports, documentation fixes, new features, and
performance improvements.

By participating in this project you agree to abide by our
[Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Please report bugs and request features via the
[GitHub issue tracker](https://github.com/waudbylab/MulticomplexNumbers.jl/issues).

When reporting a bug, please include:

- A minimal, runnable example that reproduces the problem.
- What you expected to happen and what actually happened (including any error
  messages and stack traces).
- Your Julia version (`versioninfo()`) and the package version (`]status MulticomplexNumbers`).

## Asking questions / getting support

For usage questions, please open a
[GitHub issue](https://github.com/waudbylab/MulticomplexNumbers.jl/issues) with
the `question` label. You can also consult the
[documentation](https://waudbylab.github.io/MulticomplexNumbers.jl/stable).

## Contributing code

1. **Fork** the repository and create a feature branch from `main`.
2. **Set up** the development environment:
   ```julia
   julia --project=.
   julia> using Pkg; Pkg.instantiate()
   ```
3. **Make your change.** Please follow the conventions already used in the
   codebase:
   - Keep operations type-stable and allocation-free where possible (the type is
     built on `StaticArrays`).
   - Add docstrings with an example for any new exported function.
   - Match the existing file organisation (`base.jl`, `arithmetic.jl`,
     `representations.jl`, `io.jl`).
4. **Add tests.** Every change should be covered by tests in `test/`, using
   `@testset`/`@safetestset` and `@inferred` for type-stability checks. Use
   `≈` (`isapprox`) for floating-point comparisons.
5. **Run the test suite** locally and make sure it passes:
   ```julia
   julia> using Pkg; Pkg.test()
   ```
6. **Update documentation** in `docs/src/` if you change or add user-facing
   behaviour.
7. **Open a pull request** describing the change and the motivation. CI must pass
   (tests run on Julia 1.10–1.12 across Linux, macOS, and Windows).

## Development tips

- Build the documentation locally with:
  ```bash
  julia --project=docs docs/make.jl
  ```
- Check type stability with `@code_warntype` and benchmark performance-sensitive
  code with `BenchmarkTools.jl` (see `benchmark/`).

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers this project.
