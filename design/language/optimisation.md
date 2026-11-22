# Expression optimisation

## Numerical contract

An optimisation must preserve the mathematical result and the inputs on
which it is defined. For example, `x/x` can become one if `x` is known to
be nonzero. If zero is an admissible input, keep the division.

An optimiser can establish a condition from a known value, a type or a
guarantee supplied by the caller. Such guarantees restrict the admissible
inputs and are optional. The language leaves open how to supply or check
them. Reordering floating-point operations may change their rounding.

Optimisation is optional. It may run at several stages, including after a
cell and tabulation library supplies values and after discretisation.
Reaching core normal form is a separate operation.

The defined domain means admissible input values, not the expression's spatial
domain. The examples below are grouped by the node classes involved; they are
not required optimisation passes. **FP** marks rewrites that can change
floating-point results even for finite inputs satisfying the stated conditions.
The examples assume IEEE 754 arithmetic, with real scalars unless stated otherwise.

## `Sum` and `Product`

Folding two integer literals can be exact. Regrouping a sum or factoring out a common term can change the order of rounding.

| Rewrite | Permitted? | Floating-point effect |
| ------- | ---------- | --------------------- |
| $(x + y) + z \to x + (y + z)$ | Yes | **FP**: rounding or intermediate overflow |
| $xy + xz \to x(y + z)$ | Yes | **FP**: rounding or intermediate overflow |
| `Integer(2) + Integer(3)` → `Integer(5)` | Yes | None: exact constant folding |

## `Division`

Check every denominator that a rewrite would remove. A factor of zero in front of a reciprocal still requires the reciprocal to exist.

| Rewrite | Permitted? | Floating-point effect |
| ------- | ---------- | --------------------- |
| $x/x \to 1$ | Only if $x \ne 0$ | None for finite nonzero real inputs |
| $0(1/x) \to 0$ | Only if $x \ne 0$; multiplying by zero does not remove undefinedness | **FP**: reciprocal overflow or signed zero |
| $1/3 \to 0.333$ | No: an approximation changes the mathematical value | Not applicable: forbidden by the contract |

## `Power`

The literal identities follow the defined value of a power. Combining powers with different exponents also needs a condition on the base.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $x^0 \to 1$ | `x` defined; includes $x=0$ by the core's power convention | None for finite real inputs |
| $x^1 \to x$ | Preserve the power's scalar type | None for finite real inputs |
| $x^a x^b \to x^{a+b}$ | $x>0$, real $a,b$ | **FP**: rounding, overflow or underflow |
| $(x^a)^b \to x^{ab}$ | $x>0$, real $a,b$; not a general principal-complex-power identity | **FP**: rounding, overflow or underflow |

## `Absolute` and `SquareRoot`

For real $x$, the square root of $x^2$ is its magnitude. At $x=-1$, this gives $1$.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\lvert\lvert z\rvert\rvert \to \lvert z\rvert$ | `z` may be complex | None |
| $\lvert x\rvert \to x$ | $x>0$ | None |
| $\sqrt{x^2} \to \lvert x\rvert$ | `x` real | **FP**: rounding, intermediate overflow or underflow |
| $\sqrt{x^2} \to x$ | $x\ge0$; at $x=-1$ the left side is $1$ | **FP**: rounding, intermediate overflow or underflow, or signed zero |

## `Logarithm`

A logarithm of one can be folded directly. Combining two logarithms also requires checking their domains.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\log(1) \to 0$ | Preserve the logarithm's scalar type | None for real input |
| $\log(x)+\log(y) \to \log(xy)$ | $x>0$, $y>0$; not a general principal-complex-logarithm identity | **FP**: rounding, intermediate overflow or underflow |

## `Indexed`, `Tensor`, `ComponentTensor`, `BroadcastTo` and `Zero`

A known component selection can be performed symbolically. For example, selecting entry zero of `tensor((a, b))` yields `a` when both entries have scalar value shape.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| Select an entry from an explicit `Tensor` → that entry | Valid component index | None |
| `component_tensor(a[i], (i,))[j]` → `a[j]` | Matching extents; capture-avoiding substitution | None |
| Select an entry from a broadcast scalar → that scalar | Valid component index | None |
| Select an entry from a `Zero` tensor → scalar `Zero` | Preserve the result's static type | None |

## `IndexSum`

An empty sum has no terms to evaluate. This matters when factoring out an expression that could be undefined.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\sum_i c a_i \to c\sum_i a_i$ | `c` independent of `i`; nonempty sum | **FP**: rounding or intermediate overflow |
| Sum over an index of extent zero → `Zero` | Preserve the result's static type; do not evaluate the summand | None: the specified value of an empty sum |

## `Conjugate`, `RealPart` and `ImaginaryPart`

Conjugating a value twice recovers it. When extracting a real or imaginary part, preserve the result's static type.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\overline{\overline z} \to z$ | `z` may be complex | None |
| `real_part(x)` → `x` | `x` real | None |
| `imaginary_part(x)` → `Zero` | `x` real; preserve the result's static type | None |

## `Inverse` and `Determinant`

For an invertible square matrix $A$, taking its inverse twice gives back $A$.
At a singular matrix, even the first inverse is undefined.

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $(A^{-1})^{-1} \to A$ | Square, nonsingular `A` | **FP**: numerical inversion need not recover `A` exactly |
| $\det(AB) \to \det(A)\det(B)$ | Square matrices of the same size; $AB$ denotes matrix multiplication | **FP**: rounding, overflow or underflow |

## User control and shared requirements

Users must be able to disable transformations that can affect floating-point
evaluation, preserving arithmetic operations, explicit grouping and reduction
order in UFLx expressions. Downstream consumers specify their own guarantees.

Assumptions must not be introduced silently. The [shared requirements](rewriting.md#shared-requirements) for binding,
type compatibility and [mathematical admissibility](rules.md#mathematical-admissibility)
still apply. No optimal or unique result is required.
