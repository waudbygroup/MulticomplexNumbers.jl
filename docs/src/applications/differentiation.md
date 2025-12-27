# [Numerical Differentiation](@id differentiation-applications)

```@meta
CurrentModule = MulticomplexNumbers
```

This tutorial demonstrates how to use the multicomplex step method for high-precision numerical differentiation.

---

## The Multicomplex Step Method

The multicomplex step method extends the complex step derivative to compute higher-order derivatives with machine precision accuracy, avoiding the subtractive cancellation errors that plague finite difference methods.

---

## Background: The Problem with Finite Differences

Traditional numerical differentiation uses finite differences:

```math
f'(x) \approx \frac{f(x+h) - f(x)}{h}
```

This suffers from a fundamental trade-off:
- **Large h**: Truncation error dominates (formula is only approximate)
- **Small h**: Roundoff error dominates (subtractive cancellation in `f(x+h) - f(x)`)

The optimal `h` is typically around `√ε ≈ 10⁻⁸` for Float64, limiting accuracy to about 8 digits.

---

## The Complex Step Method

For an analytic function `f`, we have the Taylor expansion:

```math
f(x + ih) = f(x) + ih \cdot f'(x) - \frac{h^2}{2}f''(x) - \frac{ih^3}{6}f'''(x) + ...
```

Taking the imaginary part:

```math
\text{Im}[f(x + ih)] = h \cdot f'(x) + O(h^3)
```

Thus:

```math
f'(x) \approx \frac{\text{Im}[f(x + ih)]}{h}
```

Since we're not computing a difference, there's **no subtractive cancellation**! We can use arbitrarily small `h`, giving nearly machine-precision derivatives.

---

## First Derivatives with `im1`

```@repl diff
using MulticomplexNumbers

# Define a test function
f(x) = sin(x) * exp(-x^2)

# Compute derivative at x = 0.5
x = 0.5
h = 1e-100  # Extremely small step - impossible with finite differences!

# Evaluate f at x + h*im1
result = f(x + h*im1)

# Extract derivative from imaginary part
f_prime = imag(result) / h

# Compute exact derivative: f'(x) = cos(x)*exp(-x²) - 2x*sin(x)*exp(-x²)
exact = cos(x)*exp(-x^2) - 2x*sin(x)*exp(-x^2)

# Error
abs(f_prime - exact)
```

The error is at the level of machine epsilon!

---

## Second Derivatives with `im1` and `im2`

For second derivatives, we use bicomplex numbers (`im1` and `im2`):

```math
f''(x) \approx \frac{\text{component}_{i_1 i_2}[f(x + h \cdot i_1 + h \cdot i_2)]}{h^2}
```

```@repl diff
using MulticomplexNumbers

# Compute second derivative of f(x) = x^4 at x = 2
# Exact answer: f''(x) = 12x² = 48

x = 2.0
h = 1e-50

# Evaluate at x + h*im1 + h*im2
result = (x + h*im1 + h*im2)^4

# The i₁i₂ component is at index 4: [1, i₁, i₂, i₁i₂]
second_derivative = component(result, 4) / h^2

# Exact value
exact = 12 * x^2

abs(second_derivative - exact)
```

---

## Third Derivatives with `im1`, `im2`, and `im3`

```@repl diff
using MulticomplexNumbers

# Compute third derivative of f(x) = x^5 at x = 1
# Exact answer: f'''(x) = 60x² = 60

x = 1.0
h = 1e-30

# Evaluate at x + h*im1 + h*im2 + h*im3
result = (x + h*im1 + h*im2 + h*im3)^5

# The i₁i₂i₃ component is at index 8
third_derivative = component(result, 8) / h^3

exact = 60.0
abs(third_derivative - exact)
```

---

## Practical Derivative Helper Function

Here's a convenient wrapper function:

```julia
using MulticomplexNumbers

"""
    derivative(f, x, n=1; h=1e-50)

Compute the n-th derivative of f at x using multicomplex step method.

# Arguments
- `f`: Function to differentiate
- `x`: Point at which to evaluate the derivative
- `n`: Order of derivative (1-4 supported)
- `h`: Step size (default: 1e-50)

# Returns
The n-th derivative of f at x

# Examples
```julia
f(x) = sin(x)
derivative(f, 1.0, 1)  # First derivative at x=1
derivative(f, 1.0, 2)  # Second derivative at x=1
```
"""
function derivative(f, x::Real, n::Int=1; h::Real=1e-50)
    if n < 1
        return f(x)
    elseif n == 1
        result = f(x + h*im1)
        return imag(result) / h
    elseif n == 2
        result = f(x + h*im1 + h*im2)
        return component(result, 4) / h^2
    elseif n == 3
        result = f(x + h*im1 + h*im2 + h*im3)
        return component(result, 8) / h^3
    elseif n == 4
        result = f(x + h*im1 + h*im2 + h*im3 + h*im4)
        return component(result, 16) / h^4
    else
        error("Derivatives of order > 4 not implemented in this example")
    end
end

# Test it
f(x) = sin(x)
x = 1.0

derivative(f, x, 1)  # f'(1) = cos(1) ≈ 0.5403
derivative(f, x, 2)  # f''(1) = -sin(1) ≈ -0.8415
derivative(f, x, 3)  # f'''(1) = -cos(1) ≈ -0.5403
derivative(f, x, 4)  # f''''(1) = sin(1) ≈ 0.8415
```

---

## Applications

### Numerical Optimization

Compute exact gradients and Hessians:

```julia
using MulticomplexNumbers

# Rosenbrock function
rosenbrock(x, y) = (1 - x)^2 + 100 * (y - x^2)^2

# Gradient at (0.5, 0.5)
x, y = 0.5, 0.5
h = 1e-50

# Partial derivatives
∂f_∂x = imag(rosenbrock(x + h*im1, y)) / h
∂f_∂y = imag(rosenbrock(x, y + h*im1)) / h

println("Gradient: ($(∂f_∂x), $(∂f_∂y))")

# Second partial: ∂²f/∂x²
∂²f_∂x² = component(rosenbrock(x + h*im1 + h*im2, y), 4) / h^2
```

### Sensitivity Analysis

Compute how outputs change with respect to parameters:

```julia
using MulticomplexNumbers

# Model: exponential decay
function decay_model(t, A, k)
    return A * exp(-k * t)
end

# Sensitivity of output to amplitude A at t=2.0
t = 2.0
A = 1.0
k = 0.5
h = 1e-50

sensitivity_A = imag(decay_model(t, A + h*im1, k)) / h
sensitivity_k = imag(decay_model(t, A, k + h*im1)) / h

println("Sensitivity to A: $(sensitivity_A)")
println("Sensitivity to k: $(sensitivity_k)")
```

### Root Finding (Newton's Method)

Newton's method with exact derivatives:

```julia
using MulticomplexNumbers

"""Newton's method using multicomplex derivatives."""
function newton(f, x0; tol=1e-10, maxiter=100)
    x = x0
    h = 1e-100

    for i in 1:maxiter
        fx = f(x)
        fpx = imag(f(x + h*im1)) / h

        if abs(fx) < tol
            return x
        end

        x = x - fx / fpx
    end

    error("Newton's method did not converge")
end

# Find root of f(x) = x^3 - 2x - 5
f(x) = x^3 - 2x - 5
root = newton(f, 2.0)
println("Root: $(root)")
println("Verification f(root) = $(f(root))")
```

---

## Practical Tips

1. **Step size**: Use `h ≈ 10⁻⁵⁰/n` for the n-th derivative. Smaller is generally better (but not required for convergence).

2. **Component indices**: For an n-th derivative, extract component `2^n` from an order-n multicomplex number:
   - 1st derivative: component 2 (index 2) from order 1
   - 2nd derivative: component 4 (index 4) from order 2
   - 3rd derivative: component 8 (index 8) from order 3
   - 4th derivative: component 16 (index 16) from order 4

3. **Computational cost**: Storage and computation scale as `2^n` for the n-th derivative. For orders above ~10, this becomes expensive.

4. **Function requirements**: Your function must:
   - Accept multicomplex arguments
   - Use operations that work on them (`+`, `-`, `*`, `/`, `^`, `exp`, `log`, `sqrt`, etc.)
   - Avoid type-specific operations or branching based on exact values

5. **Comparison with AD**:
   - Multicomplex is simpler to implement (just change the input type)
   - AD (via ForwardDiff.jl) is more efficient for gradients of functions ℝⁿ → ℝ
   - Multicomplex excels for high-order derivatives of scalar functions

---

## Limitations

1. **Analytic functions only**: The function must be complex-analytic (or at least extend to multicomplex arguments)

2. **No branching**: Avoid conditionals that depend on exact values:
   ```julia
   # Bad: branching on exact value
   f(x) = x > 0 ? sqrt(x) : 0

   # Good: smooth approximation
   f(x) = x > 0 ? sqrt(x) : sqrt(x^2 + 1e-10) - 1e-5
   ```

3. **Computational cost**: Each order doubles the component count and computational cost

4. **Very high orders**: Beyond order ~6, consider other methods (symbolic differentiation, AD with nested dual numbers, etc.)

---

## References

- **Complex Step Approximation**: [Nick Higham's blog](https://nhigham.com/2020/10/06/what-is-the-complex-step-approximation/)
- **Original Method**: Lantoine, G., Russell, R. P., Dargent, T. (2012). "Using Multicomplex Variables for Automatic Computation of High-order Derivatives." ACM TOMS 38(3):16.
- **NIST Report**: Bell, I. H., Deiters, U. K. (2021). "Precise Numerical Differentiation of Thermodynamic Functions with Multicomplex Variables." [NIST J. Res. 126:033](https://nvlpubs.nist.gov/nistpubs/jres/126/jres.126.033.pdf)

---

## See Also

- **[Creating Numbers](@ref creating)**: How to create multicomplex numbers
- **[Arithmetic Operations](@ref arithmetic)**: Mathematical operations
- **[Examples](@ref Examples)**: Quick code snippets
- **[API Reference](@ref)**: Complete function documentation
