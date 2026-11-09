# Differentiation

These reference recipes use the [shared rewrite procedure](rewriting.md#visiting-and-revisiting)
and [representation contracts](evaluation.md#representation-contracts).

## Applying derivatives

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` and applications created by its evaluation |
| Input contract         | Core expressions; applications expanded |
| Output contract        | Applications and derivatives expanded |

Prepare each derivative body under D21–D23 before selecting algebraic rules.
Substitution and dependence follow SB18, D15, SH15 and SB20–SB22 throughout.
Run `expand_applications` on exposed applications, then continue derivative
rewriting under the common procedure.

| Selection | Rules |
| --------- | ----- |
| Geometry preparation and inner derivatives | D21–D23, D13 |
| Independent expression | L10 |
| Variable and application derivatives | D11–D12 |
| Derivative fallback | D01–D02, D19–D20, D03–D10, L04, L01, ST01, L02–L03, in order, with $\delta b = \mathrm{D} b$ |

### Geometry preparation

A variable may occur in a map’s parameterisation even when absent from the
mapped expression’s value operands. Check construction dependencies before
selecting constant or commuting rules.
Prepare the derivative's body, preserving its type under SH02–SH04; its
variable declaration, point and direction are not preparation targets.

| Rule | Condition | Requirement |
| ---- | --------- | ----------- |
| D21 | Inner derivatives remain | Resolve them inside out, each with its own variable, point and direction |
| D22 | A variable occurs in a non-value dependency | Expose the dependence through integral expansion, spatial mapping and library reconstruction, or use an explicit derivative rule. Recheck exposed expressions; independent geometry needs no expansion. |
| D23 | Dependence remains unresolved or no compatible fixed-space representation is supplied | Raise `ValueError` identifying the dependency. Neither infer independence from value operands alone nor choose implicit transport. |

Preparation uses the spatial recipes, resolving newly exposed
derivatives by the same procedure. Prepared-geometry and progress requirements
apply; local algebraic rules do not invoke whole passes. Library derivative
rules must retain dependencies needed by later derivatives or reject them.

### Directional derivative rules

| Rule | Expression $f$ | Rewrite of the derivative at $a$ |
| --- | --- | --- |
| D11 | $f=x$ for a single function, or $f=x_k$ for a tuple variable | $h$ for a single function, or $h_k$ for the matching tuple member |
| D12 | `Apply(b, G)` | `apply(Db, G_at)` plus one `apply(b_at, G_at_with_Dg)` per slot: replace just that slot of `G_at` by its derivative. Differentiate tuple-valued slots member by member. Here `b_at` and `G_at` are evaluated at $a$; inserted directions remain unchanged. |
| D13 | `DirectionalDerivative` | Resolve the inner derivative, including D21–D23, before differentiating its result. Its point and direction are outside its variable scope, so their free dependencies participate in the outer derivative. |

| Rule | Applies to | Requirement |
| --- | --- | --- |
| D15 | Dependence on the derivative variable | First use the construction walk to detect dependencies, including geometry and spaces; apply D21–D23 before SH15 and SB22. Tuple variables vary simultaneously, including under conjugation. Dual points and directions may be symbolic forms whose applications the recipe evaluates. |

For real scalar functions `u` and `a` in a space of constants,
`directional_derivative(u * u, u, a, u)` rewrites to `u * a + a * u`.

For a real scalar function `u` in a space of constants and a constant array
`a` of shape `(m,)`, `directional_derivative(u, u, a[i], Integer(1))` lowers
to `indexed(broadcast_to(Integer(1), (m,)), (i,))`. Its value is one for every `i`.
Type completion retains indices even from unused members of tuple points
and directions.

| Rule | Expression $f$ | Rewrite of the derivative at $a$ |
| --- | --- | --- |
| L10 | $f$ passes the directional independence check SH15 for $x$ | Zero with the type of the derivative |

| Rule | Applies to | Requirement |
| --- | --- | --- |
| SB18 | Derivative substitution | After D21–D23, for `DirectionalDerivative(b, x, a, h)`, match original operands before substituting `[x := a]` into undifferentiated value operands; insert `h` unchanged. Remaining domains, maps and spaces must be independent of $x$ or handled by an explicit derivative rule. |

## Applying gradients

Expand compound expressions into component derivatives, retaining terminal
and mapped gradients. For a scalar terminal field `u` on `P`,

$$
\partial_j(u u) = (\partial_j u)u + u(\partial_j u).
$$

Each $\partial_j u$ comes from the retained `gradient(u, P)`; binding `j`
reconstructs the result. A terminal gradient alone remains whole.

| `apply_gradients(f)` |  |
| -------------------- | - |
| Evaluates            | `Gradient`, `Divergence` and `Curl` |
| Input contract       | Core expressions; applications and derivatives expanded |
| Output contract      | Applications, derivatives and gradients expanded |

After spatial rewriting, run `expand_applications` on exposed applications
and continue under the common procedure. Gradients exposed by a later
pull-back expansion belong to the next call to this recipe; `lower` schedules
that call in its next spatial round.

### Whole spatial operations

Choose among these cases before requesting any component derivatives.

| Selection, in order | Rules |
| --- | --- |
| Move an application outside the gradient | ST07 |
| Keep Piola patterns and retained gradients | CC06–CC07, ST12 |
| Expand divergence and curl | D16 |
| Reconstruct the gradient from components | ST08 |

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| ST07 | `Gradient(Apply(B, G), P)` | `Apply(Gradient(B, P), G)` |
| D16 | `Divergence` or `Curl` not matching CC06 or CC07 | Its composition of gradients from [Divergence and Curl](geometry.md#divergence-and-curl) |
| ST08 | `Gradient(f, P)`, $f$ not an `Apply` and the gradient does not match ST12 | `component_tensor(∂_j f[A], A + (j,))`, using the component rules selected by [Applying gradients](#applying-gradients) |

Gradient reconstruction uses derived indices `A` and `j`, with one index in
`A` per value axis of the operand. Empty lists use the factory identities.

### Retained gradients

ST12 defines the retained gradients admitted by the gradient contract.

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| ST12 | `Gradient(b, P)`, $b$ with a domain and either a terminal, `Gradient`, `PushForward`, node matching CC06 or CC07, or `PullBack` by a cell or `SideMap`; exclude `Coordinate(P)` when $P$ has its ambient dimension | Retain the whole gradient unchanged |

### Retained Piola patterns

These divergence and curl patterns take precedence over generic expansion.
Here `k` and `n` are the topological and ambient dimensions of the mesh `D`.
Their cell pull-backs use [Piola identities](lowering.md#piola-identities).

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| CC06 | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)` | Remains for the Piola rule of a cell pull-back |
| CC07 | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$ | Remains for the Piola rule of a cell pull-back |

### Component derivatives

Resolve $\partial_j b$ recursively using these cases; ST09 stops at a retained
gradient. The request has `b`’s shape and an additional free index `j`.

| Selection, in order | Rules |
| --- | --- |
| Constants and coordinates | L11, D17 |
| Components of retained gradients | ST09 |
| Chain rule and applications | CC01, ST10 |
| Shared derivative fallback | D01–D02, D19–D20, D03–D10, L04, L01, ST01, L02–L03, in order, with $\delta b = \partial_j b$ |

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| L11 | $\partial_j b$, $b$ passes the spatial independence check SH16 | Zero with $b$'s shape and argument spaces, and its free indices followed by $j$ |
| D17 | $\partial_j$ `Coordinate(P)`, $P$ of its ambient dimension | `component_tensor(I[i, j], (i,))`, with $I$ the identity of shape `(n, n)` and a derived $i$ |
| ST09 | $\partial_j b$, with `Gradient(b, P)` matching ST12 | `component_tensor(gradient(b, P)[B + (j,)], B)`, with one derived index in `B` per value axis of $b$; ST12 retains the introduced gradient |
| CC01 | $\partial_j$ `PullBack(g, φ)`, φ a `DomainMap` with parameter domain $P$ | `component_tensor(index_sum(H[A + (l,)] * J[l, j], (l,)), A)`, with `H = pull_back(gradient(g, Q), φ)`, `J = jacobian(φ)`, and $Q$ the ambient domain of φ |
| ST10 | $\partial_j$ `Apply(B, G)` | `Apply(∂_j B, G)`, with every supplied value unchanged |

In the chain rule, `A` indexes the value axes of `g` and `l` is a distinct
contraction index; scalar `g` uses `A == ()`. On embedded meshes the gradient
is tangential. The coordinate gradient remains a projection until the
pull-back to a cell.

| Rule | Applies to | Requirement |
| --- | --- | --- |
| SH24 | A replacement contains $\partial_j b$ | Resolve this component-derivative request with the component rows of `apply_gradients`. It is notation, not a core node. ST08 creates these requests; ST09 ends a request at a retained gradient. No unresolved request may appear in the returned expression. |

## Derivatives of the node classes

Directional and spatial derivatives share these rules. Their specific rules
take precedence. The linear derivative rules are in [Linearity](rewriting.md#linearity).
Extensions can add derivative rules to either group.

| Rule | Expression | Derivative |
| --- | --- | --- |
| D01 | `Product(b, c)` | $\delta b \, c + b \, \delta c$ |
| D02 | `Division(b, c)` | $\delta b / c - (b \, \delta c) / (c c)$ |
| D19 | `Power(b, c)`, $c$ an `Integer`, `Float` or `Complex` literal with value zero, or a scalar `Zero` | Zero with the type of the derivative |
| D20 | `Power(b, c)`, $c$ an `Integer`, `Float` or `Complex` literal with value one | $\delta b$, with type completion against the derivative |
| D03 | `Power(b, c)`, independence of $c$ established by SH15 or SH16, or an analysis under SH18; exclude D19 and D20 | $c \, b^{c - 1} \, \delta b$ |
| D04 | `Power(b, c)`, independence of $c$ not established | $b^c (\log(b) \, \delta c + c \, \delta b / b)$ |
| D05 | `Absolute(b)` | $\operatorname{Re}(\bar b \, \delta b) / \lvert b \rvert$, with `real_part`, so the derivative is real |
| D06 | `SquareRoot(b)` | $\delta b / (2 \sqrt b)$ |
| D07 | `Logarithm(b)` | $\delta b / b$ |
| D08 | `Determinant(A)`, shape `(n, n)` | $\sum_{i,j=0}^{n-1} C_{ij}(A)\,\delta A_{ij}$, with cofactors as defined below |
| D09 | `Inverse(A)`, equal axis extents, and $A$ invertible | `-inverse(A) @ delta_A @ inverse(A)`, with `delta_A` = $\delta A$ |
| D10 | `Inverse(A)`, $A$ real with shape `(n, k)`, $k<n$, and full column rank | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |

D08 uses $C_{ij}(A)=(-1)^{i+j}\det(A_{\widehat i,\widehat j})$, where the
minor removes row $i$ and column $j$. For `n >= 2`, construct each minor
with `Tensor` and integer component selections, and expand the finite sum
with scalar addition and multiplication. For `n == 1`, the result is
`delta_A[0, 0]`. For `n == 0`, it is zero with the type of the derivative.

Literal values in D19 and D20 are compared exactly, including signed zero;
`Complex` has that value when its imaginary part is zero. These cases apply
before substituting the differentiation point. A variable exponent that happens
to be zero or one at that point still uses D04. An undefined base remains
undefined under mathematical admissibility.

The constant-exponent power rule avoids a logarithm, permitting negative real
bases where differentiable. Results are core expressions.

## Nested derivative examples

Displayed results omit zero terms and use schematic fresh names;
rewriting need not simplify algebra.

Let `u`, `a`, `b`, `h` and `k` be real scalar functions in a constant space,
with `a`, `b`, `h` and `k` independent of `u`.

| Input | Rule trace and result |
| ----- | --------------------- |
| `directional_derivative(directional_derivative(u * u, u, u, h), u, a, k)` | The inner derivative gives `h * u + u * h` by D01, D11 and SB18. Its point `u` is free under SB04. The outer derivative gives `h * k + k * h`. |
| `directional_derivative(directional_derivative(u * u, u, b, h), u, a, k)` | The inner derivative gives `h * b + b * h`. SH15 then finds no free `u`, so L10 gives the typed zero. |
