# Substitution

## Substitution and binding

Substitution preserves bound names and keeps inserted free names free.
Here `j` is fresh:

| Substitution | Result |
| ------------ | ------ |
| `replace(index_sum(a[i] * c, (i,)), {c: b[i]})` | `index_sum(a[j] * b[i], (j,))`; `i` remains free |
| `replace(form(a, (v,)), {v: w})` | The original form; `v` is bound |
| `replace(DirectionalDerivative(f, u, u, h), {u: u0})` | `DirectionalDerivative(f, u, u0, h)`; the point is outside the binding scope |

### Free occurrences

A name is an `Index` or a `Function`. Free occurrence uses the same binding
rules for two walks, chosen by the operation.

| Rule | Walk or condition | Definition |
| ---- | ----------------- | ---------- |
| SB20 | Free occurrence of a name | Encounter a structurally equal name on the chosen walk with no enclosing binder binding it in that operand. Binding declarations do not themselves count as free occurrences. Names mentioned only in labels do not count. |
| SB21 | Construction walk, used by substitution and capture checks | Follow construction dependencies, including domains, maps, spaces and the spaces of declared functions, respecting SB02–SB04. A declaration is not a free use of its own name; its space dependencies remain free. Label reconstruction follows [Derived indices](expressions.md#derived-indices), but label contents add no free occurrences. |
| SB22 | Value walk, used after differentiation dependencies are exposed | Follow expression value operands SB01 and index operands, respecting SB02–SB04. A `Function` contributes itself; do not enter its space or label. Do not enter domains, maps or other non-expression dependencies. |

Tuples are walked member by member. A name may be bound in one operand and
free in another. Unless stated otherwise, substitution and capture avoidance
use the construction walk. Directional differentiation uses it to detect
hidden dependence before using the value walk on the prepared expression.

| Example | Free occurrences |
| ------- | ---------------- |
| `form(e, (u,))` | Occurrences of `u` in `e` are bound |
| `DirectionalDerivative(f, u, u, h)` | `u` is bound in `f` but free in the point operand |
| `Function(V, ("label", u))` | The label introduces no free `u` in either walk |
| `pull_back(f, φ)`, with `u` absent from `f` but free in `φ.parameterisation` | `u` occurs in the construction walk but not the value walk |

The static `free_indices` attribute retains its definition in
[Free indices](expressions.md#free-indices). The value walk is a dependency
check, not a replacement for that attribute.

### Scopes and substitution

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB01 | Value operands and labels | Value operands are `AbstractExpression` construction arguments, including those in tuples. Terminals have none. Labels identify objects and follow [Derived indices](expressions.md#derived-indices). Only these nodes bind names. |

| Rule | Binding node | Value operands | Bound names | Bound in |
| ---- | ------------ | -------------- | ----------- | -------- |
| SB02 | [`ComponentTensor(e, indices)`, `IndexSum(e, indices)`](expressions.md#componenttensor-and-indexsum) | `e` | The indices of `indices` | `e` and `indices` |
| SB03 | [`Form(e, arguments)`](forms.md#form) | `e` | The functions of `arguments` | `e` and `arguments` |
| SB04 | [`DirectionalDerivative(f, x, a, h)`](geometry.md#directionalderivative) | `f`, `a`, `h` | $x$ | $f$ and $x$ |

All binders use lexical scope. Renaming preserves their ordered signatures
(function spaces or index extents). Declarations contribute type dependencies,
not free uses of their names. Alpha-equivalence compares declarations while
preserving free identities, for both functions and indices.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB05 | Substitution traversal | Use the construction walk SB21 and the binding scopes SB02–SB04. Rebuild derived labels under [Derived indices](expressions.md#derived-indices). Each pass states which operands it transforms. |
| SB06 | Fresh bound names | Making a binding node *fresh* for expressions $s_1, \dots, s_n$ means renaming each bound name that occurs free in any $s_k$ under SB20–SB21. Rename it only within its binding scope. This prevents an inserted free name from becoming bound accidentally. |

Substitution replaces $r_1, \dots, r_n$ by $s_1, \dots, s_n$
simultaneously. Within a binding scope, follow these steps.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB07 | Within a binding scope, step 1 | Drop each pair whose key $r_k$ has a bound name free in it under SB20–SB21. |
| SB08 | Within a binding scope, step 2 | Make the binding node fresh for the remaining replacement values $s_k$. |
| SB09 | Within a binding scope, step 3 | Substitute the remaining pairs into the operands in that scope. |
| SB10 | Outside the binding scope | Operands outside the binding scope use the original pairs. |

`replace(e, mapping)` substitutes each key of the dictionary `mapping` with
its value, simultaneously.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB13 | Matching | Match whole original nodes before their dependencies, respecting binding scopes and label rules. The outermost eligible match wins, even when keys overlap. Do not search inserted values or rebuilt nodes. |
| SB14 | Replacement order | Order replacement pairs by the structural key order of [`Base`](rules.md#the-base-object-base), preserving that order through scope filtering and in fresh-name labels. Dictionary insertion order has no effect. |

```python
replace(a + b, {a + b: c, a: d})  # c
replace(a + b, {a: d, d + b: c})  # d + b; no match after rebuilding
replace(a, {a: b, b: c})  # b; inserted values are not searched
```

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB15 | Compositions | Substitution commutes with a fixed composition when the static-type choice is unchanged. For example, replacing `A` by `B` in `transpose(A)` gives `transpose(B)`, with the derived indices of `B`. |

### Constructing fresh names

To make names fresh for $s_1, \ldots, s_n$, let `body` be the tuple of operands
in the binder's scope. Use these names.

| Rule | Bound name | Fresh name |
| ---- | ---------- | ---------- |
| SB11 | Index `i` | `Index(("uflx.rename", i, body, (s_1, ..., s_n)))` |
| SB12 | Function `v` in $V$ | `Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))` |

The [derived labels](expressions.md#derived-indices) keep these names distinct
and absent from `body` and the inserted expressions.

## Capture avoidance example

Displayed results omit zero terms and use schematic fresh names;
rewriting need not simplify algebra.

Let `u` and `u0` be real scalar functions in a constant space, and `a` and `b`
be constant arrays. The index `i` has the same extent in both arrays.

| Input | Rule trace and result |
| ----- | --------------------- |
| `directional_derivative(index_sum(u * a[i], (i,)), u, u0, b[i])` | SB06–SB09 rename the sum's bound index to fresh `j` before ST01 moves differentiation inside. D01 and D11 give `index_sum(b[i] * a[j], (j,))`. The direction's `i` stays free. |

For `a = (2, 3)` and `b = (7, 11)`, the result has values `(35, 55)` over
free `i`. Capturing `i` would instead give the scalar `47`.
