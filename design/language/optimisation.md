# Expression optimisation

Expression optimisation rewrites an expression to reduce its expected evaluation
cost, for example by removing redundant calculations or factorising terms.

Expression optimisation is optional and separate from establishing core normal
form. It may be applied at multiple stages, including after discretisation and
substitution of values supplied by a cell and tabulation library.

Optimisation keeps the same mathematical result. Where the original expression
is undefined, the result must also be undefined. The order and rounding of
floating-point operations may change.

Some simplifications need extra information. For example, `x/x` can become `1`
only when `x` is known to be nonzero. This may follow from a known value, the
expression's type, or a guarantee supplied by the caller. Such guarantees
restrict the inputs for which the optimiser promises to preserve the result.
Without enough information, the expression is left unchanged. Callers need not
supply guarantees, and the language does not yet prescribe how to supply or check
them.

The defined domain means admissible input values, not the expression's spatial
domain. The examples below are grouped by the node classes involved; they are
not required optimisation passes. **FP** marks rewrites that can change
floating-point results even for finite inputs satisfying the stated conditions.
The examples assume IEEE 754 arithmetic, with real scalars unless stated otherwise.

## `Sum` and `Product`

| Rewrite | Permitted? | Floating-point effect |
| ------- | ---------- | --------------------- |
| $(x + y) + z \to x + (y + z)$ | Yes | **FP**: rounding or intermediate overflow |
| $xy + xz \to x(y + z)$ | Yes | **FP**: rounding or intermediate overflow |
| `Integer(2) + Integer(3)` → `Integer(5)` | Yes | None: exact constant folding |

## `Division`

| Rewrite | Permitted? | Floating-point effect |
| ------- | ---------- | --------------------- |
| $x/x \to 1$ | Only if $x \ne 0$ | None for finite nonzero real inputs |
| $0(1/x) \to 0$ | Only if $x \ne 0$; multiplying by zero does not remove undefinedness | **FP**: reciprocal overflow or signed zero |
| $1/3 \to 0.333$ | No: an approximation changes the mathematical value | Not applicable: forbidden by the contract |

## `Power`

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $x^0 \to 1$ | `x` defined; includes $x=0$ by the core's power convention | None for finite real inputs |
| $x^1 \to x$ | Preserve the power's scalar type | None for finite real inputs |
| $x^a x^b \to x^{a+b}$ | $x>0$, real $a,b$ | **FP**: rounding, overflow or underflow |
| $(x^a)^b \to x^{ab}$ | $x>0$, real $a,b$; not a general principal-complex-power identity | **FP**: rounding, overflow or underflow |

## `Absolute` and `SquareRoot`

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\lvert\lvert z\rvert\rvert \to \lvert z\rvert$ | `z` may be complex | None |
| $\lvert x\rvert \to x$ | $x>0$ | None |
| $\sqrt{x^2} \to \lvert x\rvert$ | `x` real | **FP**: rounding, intermediate overflow or underflow |
| $\sqrt{x^2} \to x$ | $x\ge0$; at $x=-1$ the left side is $1$ | **FP**: rounding, intermediate overflow or underflow, or signed zero |

## `Logarithm`

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\log(1) \to 0$ | Preserve the logarithm's scalar type | None for real input |
| $\log(x)+\log(y) \to \log(xy)$ | $x>0$, $y>0$; not a general principal-complex-logarithm identity | **FP**: rounding, intermediate overflow or underflow |

## `Indexed`, `Tensor`, `ComponentTensor`, `BroadcastTo` and `Zero`

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| Select an entry from an explicit `Tensor` → that entry | Valid component index | None |
| `component_tensor(a[i], (i,))[j]` → `a[j]` | Matching extents; capture-avoiding substitution | None |
| Select an entry from a broadcast scalar → that scalar | Valid component index | None |
| Select an entry from a `Zero` tensor → scalar `Zero` | Preserve the result's static type | None |

## `IndexSum`

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\sum_i c a_i \to c\sum_i a_i$ | `c` independent of `i`; nonempty sum | **FP**: rounding or intermediate overflow |
| Sum over an index of extent zero → `Zero` | Preserve the result's static type; do not evaluate the summand | None: the specified value of an empty sum |

## `Conjugate`, `RealPart` and `ImaginaryPart`

| Rewrite | Condition | Floating-point effect |
| ------- | --------- | --------------------- |
| $\overline{\overline z} \to z$ | `z` may be complex | None |
| `real_part(x)` → `x` | `x` real | None |
| `imaginary_part(x)` → `Zero` | `x` real; preserve the result's static type | None |

## `Inverse` and `Determinant`

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
