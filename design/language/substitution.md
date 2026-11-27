# Substitution

## Substitution and binding

If a substitution would insert a free `x` inside a form that binds `x`,
give the form argument a fresh name first. If nothing is inserted there,
keep the argument name.

Take distinct scalar functions `w` and `x` in the same constant space.

| Expression | Result of `replace(expression, {w: x})` |
| ---------- | --------------------------------------- |
| `form(x, (x,))` | The original form. There is no `w` to replace. |
| `form(w * x, (x,))` | `form(x * y, (y,))` for fresh `y`. The inserted `x` stays free. |

The same rule protects an index. Replacing `c` by `b[i]` in
`index_sum(a[i] * c, (i,))` gives `index_sum(a[j] * b[i], (j,))`
for fresh `j`. The sum still runs over the entries of `a`.

### Free occurrences

A name (an `Index` or a `Function`) is free at an occurrence outside the
scope of its declaration. In `DirectionalDerivative(f, u, u, h)`, the
derivative binds `u` in `f`, but the evaluation point is a free use of
`u`. If `h` is independent of `u`, replacing `u` by `u0` changes the
point and gives `DirectionalDerivative(f, u, u0, h)` for a fresh `u0`
in the same space.

Which occurrences can be seen depends on which edges of the graph are
followed. For example, `u` might appear only in the parameterisation of
the map in `pull_back(f, φ)`. The construction walk sees that dependency.
Directional differentiation must expose it before using the value walk.

| Rule | Walk or condition | Definition |
| ---- | ----------------- | ---------- |
| SB20 | Free occurrence | Look for a structurally equal name on the chosen walk. Ignore occurrences bound in that operand, the declaration of the name itself, and names mentioned only inside labels. |
| SB21 | Construction walk | Follow construction dependencies, including domains, maps, spaces and the spaces of declared functions. Respect the scopes below. Rebuild [derived labels](expressions.md#derived-indices) even though their contents add no free occurrences. |
| SB22 | Value walk | Follow expression value operands and index operands within their scopes. A `Function` contributes itself. Leave spaces, labels, domains, maps and other non-expression metadata untouched. |

The space of a declared function can depend on other free names.
Both walks visit tuples member by member. Unless otherwise stated, substitution
and capture checks use the construction walk. The static
[`free_indices` attribute](expressions.md#free-indices) has its own definition.

### Scopes and substitution

A binder declares which names it owns and where it binds them. For example,
`replace(form(v, (v,)), {v: w})` leaves the form alone because the `v`
in its body is bound.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB01 | Value operands and labels | Value operands are `AbstractExpression` construction arguments, including tuple members. Terminals have none. Labels identify objects and follow [Derived indices](expressions.md#derived-indices). Only the nodes below bind names. |

| Rule | Binding node | Value operands | Bound names | Bound in |
| ---- | ------------ | -------------- | ----------- | -------- |
| SB02 | [`ComponentTensor(e, indices)`, `IndexSum(e, indices)`](expressions.md#componenttensor-and-indexsum) | `e` | The declared indices | `e` and `indices` |
| SB03 | [`Form(e, arguments)`](forms.md#form) | `e` | The functions in `arguments` | `e` and `arguments` |
| SB04 | [`DirectionalDerivative(f, x, a, h)`](geometry.md#directionalderivative) | `f`, `a`, `h` | $x$ | $f$ and $x$ |

Renaming a bound name updates occurrences in its lexical scope and preserves
the ordered binding signature (function spaces or index extents).
Alpha-equivalence allows consistent renaming of bound functions and indices
while preserving free identities.

#### Matching the original expression

`replace(e, mapping)` makes all replacements simultaneously. Test a whole
node before looking inside it. Once it matches, insert the supplied value
and finish with that occurrence.

```python
replace(a + b, {a + b: c, a: d})  # c
replace(a + b, {a: d, d + b: c})  # d + b
replace(a, {a: b, b: c})          # b
```

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB13 | Overlapping keys | The outermost eligible match in the original graph wins. Respect scopes and label rules. Inserted values and rebuilt nodes are never searched for further matches. |
| SB14 | Order | Keep pairs in the [structural order of their keys](rules.md#the-base-object-base), including when filtering by scope or constructing fresh-name labels. Dictionary insertion order has no effect. |

#### Substituting under a binder

When replacing `w` by `x` in a sum of forms, decide separately for each term
whether an insertion would capture `x`. The example at the start of the
chapter shows the two possible outcomes.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB05 | Traversal | Use the construction walk and the binding scopes above. Rebuild derived labels. Each pass states which operands it transforms. |
| SB06 | Making a binder fresh for given values | Rename each bound name that occurs free in any of those values on the construction walk. Rename it only within its scope. |
| SB07 | Choose eligible keys | In this scope, drop any key that contains a free use of a name declared by this binder. |
| SB08 | Choose fresh names | Find the eligible keys that actually match in this scope, respecting outermost matches and nested scopes. Ignore matches inside labels for this step. Make the binder fresh for the corresponding values in structural key order. |
| SB09 | Rebuild the scope | Apply the eligible substitutions to the original scoped operands, including label reconstruction under SB05. Rename the original bound occurrences. |
| SB10 | Visit operands outside the scope | Use the original mapping, including pairs dropped for the scoped operands. |

A key with no eligible match has no effect. For a fixed composition such as
`transpose(A)`, replacement also updates indices derived from the original
operands.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB15 | Compositions | When the static-type choice stays the same, substitution commutes with the composition. Replacing `A` by `B` in `transpose(A)` gives `transpose(B)`, including its derived indices. |

### Constructing fresh names

The `y` in the opening example stands for a reproducible function in the
same space as `x`. Its label records the old name, the tuple `body` of
operands in the binder's scope and the ordered values $s_1,\ldots,s_n$
for which the name must be fresh.

| Rule | Bound name | Fresh name |
| ---- | ---------- | ---------- |
| SB11 | Index `i` | `Index(("uflx.rename", i, body, (s_1, ..., s_n)))` |
| SB12 | Function `v` in $V$ | `Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))` |

The [derived labels](expressions.md#derived-indices) keep the fresh names
distinct and absent from `body` and the inserted expressions.

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
