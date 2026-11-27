# Differentiation

## Applying derivatives

Differentiate the expression at the supplied point in the supplied direction.
For a real scalar function `u`, the product rule gives

```python
directional_derivative(u * u, u, a, h)  # h * a + a * h
```

Both undifferentiated copies of `u` take the value `a`. The direction
`h` enters as supplied.

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` and applications created by its evaluation |
| Input contract         | Core expressions; applications expanded |
| Output contract        | Applications and derivatives expanded |

First resolve inner derivatives and expose any dependence on geometry. Then
check for independence and apply the algebraic rules in the order below.
Expand any applications they expose and continue rewriting. The
[directional derivative rules](#directional-derivative-rules) specify
which operands take the value at the differentiation point.

| Selection | Rules |
| --------- | ----- |
| Geometry preparation and inner derivatives | D21–D23, D13 |
| Independent expression | L10 |
| Variable and application derivatives | D11–D12 |
| Derivative fallback | B02, D01–D02, D19–D20, D03–D10, L04, L01, ST01, L02–L03, in order, with $\delta b = \mathrm{D} b$ |

### Geometry preparation

An integral of the constant field `Integer(1)` can depend on the variable
being differentiated. If the coordinates of mesh `D` depend on `u`,
then `integral(Integer(1), D)` measures a length, area or volume that
can change with `u`. Check the stored geometry before declaring the
derivative zero.

Prepare only the derivative body, preserving its [type](rewriting.md#type-compatibility-of-rewrites).
Its variable declaration, point and direction stay as supplied.

| Rule | Condition | Requirement |
| ---- | --------- | ----------- |
| D21 | Inner derivatives remain | Resolve them inside out, each with its own variable, point and direction |
| D22 | A variable occurs in a non-value dependency | Expose the dependence through integral expansion, spatial mapping and library reconstruction, or use an explicit derivative rule. Recheck exposed expressions; independent geometry needs no expansion. |
| D23 | Dependence remains unresolved or no compatible fixed-space representation is supplied | Raise `ValueError` identifying the dependency. Neither infer independence from value operands alone nor choose implicit transport. |

Preparation uses the spatial recipes, resolving newly exposed
derivatives by the same procedure. Prepared-geometry and progress requirements
apply; local algebraic rules do not invoke whole passes. Library derivative
rules must retain dependencies needed by later derivatives or reject them.

Geometry-dependent forms require an explicit identification of their
argument spaces across mesh motion. Preparation preserves the input
signature and output domain. Mapped spaces identify their argument DOFs through the fixed
reference space. A coordinate variation changes the physical basis.
CC10 exposes the mapping and CC11 reads the independent reference data.
For a space with physical nodal data and geometry-dependent basis
reconstruction, hold the supplied nodal data fixed and differentiate the
reconstruction.

Libraries must expose geometry-dependent reconstruction, including basis and
coefficient transformations, as symbolic operands or derivative rules.

Independent field data, topology, connectivity, index extents and discrete
element choices remain fixed. Geometry must remain admissible and
differentiable; topology changes and singular or nonsmooth configurations
are excluded.

### Directional derivative rules

Only occurrences of the variable in the body are bound. The direction can
even be `u` itself. For real scalars in a constant space,
`directional_derivative(u * u, u, a, u)` becomes `u * a + a * u`.

| Rule | Expression $f$ | Rewrite of the derivative at $a$ |
| --- | --- | --- |
| D11 | $f=x$ for a single function, or $f=x_k$ for a tuple variable | $h$ for a single function, or $h_k$ for the matching tuple member |
| D12 | `Apply(b, G)` | `apply(Db, G_at)` plus one `apply(b_at, G_at_with_Dg)` per slot: replace just that slot of `G_at` by its derivative. Differentiate tuple-valued slots member by member. Here `b_at` and `G_at` are evaluated at $a$; inserted directions remain unchanged. |
| D13 | `DirectionalDerivative` | Resolve the inner derivative, including D21–D23, before differentiating its result. Its point and direction are outside its variable scope, so their free dependencies participate in the outer derivative. |

| Rule | Applies to | Requirement |
| --- | --- | --- |
| D15 | Dependence on the derivative variable | First use the construction walk to detect dependencies, including geometry and spaces; apply D21–D23 before SH15 and SB22. Tuple variables vary simultaneously, including under conjugation. Map-valued points and directions may be symbolic forms with admissible signatures whose applications the recipe evaluates. |

For a fixed bilinear map $B$ and fixed $v$, varying $u$ in direction $h$
in $B(u,v)$ gives $B(h,v)$. If both slots vary, D12 adds one term for
each slot.

Fixed-direction Hessians use independent direction functions. Mixed
derivatives commute only under the required smoothness and fixed-direction
hypotheses.

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

Process the whole operator first. For example,
`gradient(apply(A, G), P)` becomes `apply(gradient(A, P), G)`.
A cell pull-back can also use a whole divergence or curl directly.
Expand into components only after checking these cases.

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

Leave the gradient whole when the remaining operation belongs to the
field's library or needs a map. For example, `gradient(u, P)` stays as
written for a terminal function `u` on `P`.

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| ST12 | `Gradient(b, P)`, $b$ with a domain and either a terminal, `Gradient`, `PushForward`, node matching CC06 or CC07, or `PullBack` by a cell or `SideMap`; exclude `Coordinate(P)` when $P$ has its ambient dimension | Retain the whole gradient unchanged |

### Retained Piola patterns

A mapped field can carry a divergence or curl by a Piola identity. For
example, retain the divergence of a contravariant Piola field so a cell
pull-back can use its reference divergence directly. These patterns take
precedence over generic expansion.
Here `k` and `n` are the topological and ambient dimensions of the mesh `D`.
Their cell pull-backs use [Piola identities](lowering.md#piola-identities).

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| CC06 | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)` | Remains for the Piola rule of a cell pull-back |
| CC07 | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$ | Remains for the Piola rule of a cell pull-back |

### Component derivatives

Compute a spatial derivative in direction `j` and keep `j` free. In the
full-dimensional coordinate field `X = coordinate(P)`, the derivative of
`X[i]` in direction `j` is the identity entry $I_{ij}$. For a scalar
terminal field `u`, the derivative in direction `j` is `gradient(u, P)[j]`.
The component request has the shape of its operand and one additional
free index `j`.

| Selection, in order | Rules |
| --- | --- |
| Constants and coordinates | L11, D17 |
| Components of retained gradients | ST09 |
| Chain rule and applications | CC01, ST10 |
| Shared derivative fallback | B02, D01–D02, D19–D20, D03–D10, L04, L01, ST01, L02–L03, in order, with $\delta b = \partial_j b$ |

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

Directional and spatial differentiation use the same algebraic identities
after their specific rules have been tried. The symbol $\delta b$ stands
for the derivative of operand $b$. In a directional derivative, evaluate
undifferentiated operands at the supplied point. The
[linearity rules](rewriting.md#linearity) cover linear operations.
Extensions can add derivative rules to either group.

### Products, quotients and conditionals

A numerical conditional takes the derivative of the selected branch. For
example, differentiating `conditional(p, u, v)` selects between the
derivatives of `u` and `v`. Boolean values have no directional or spatial
derivative, and the predicate itself stays undifferentiated.

| Rule | Expression | Derivative |
| --- | --- | --- |
| B02 | `Conditional(p, a, b)`, numerical | `conditional(p_at, delta_a, delta_b)`; for a directional derivative, `p_at` is `p` with the differentiation point substituted under SB18; for a spatial derivative it is `p` unchanged |
| D01 | `Product(b, c)` | $\delta b \, c + b \, \delta c$ |
| D02 | `Division(b, c)` | $\delta b / c - (b \, \delta c) / (c c)$ |

### Powers

For a fixed exponent, use $c b^{c-1}\,\delta b$ wherever that formula is
defined. A literal zero exponent gives one at every defined base, including
zero, so its derivative is zero there. For a symbolic exponent whose
independence has been established, the same formula assumes
[mathematical admissibility](rules.md#mathematical-admissibility). If the
base and this exponent both evaluate to zero, it contains $0^{-1}$.

| Rule | Expression | Derivative |
| --- | --- | --- |
| D19 | `Power(b, c)`, $c$ an `Integer`, `Float` or `Complex` literal with value zero, or a scalar `Zero` | Zero with the type of the derivative |
| D20 | `Power(b, c)`, $c$ an `Integer`, `Float` or `Complex` literal with value one | $\delta b$, with type completion against the derivative |
| D03 | `Power(b, c)`, independence of $c$ established by SH15 or SH16, or an analysis under SH18; exclude D19 and D20 | $c \, b^{c - 1} \, \delta b$ |
| D04 | `Power(b, c)`, independence of $c$ not established | $b^c (\log(b) \, \delta c + c \, \delta b / b)$ |

Literal values in D19 and D20 are compared exactly, including signed zero;
`Complex` has that value when its imaginary part is zero. These cases apply
before substituting the differentiation point. A variable exponent that happens
to be zero or one at that point still uses D04. An undefined base remains
undefined under mathematical admissibility.

The constant-exponent power rule avoids a logarithm, permitting negative real
bases where differentiable. Results are core expressions.

### Absolute value, square root and logarithm

For real $b>0$, differentiating $\sqrt b$ gives $\delta b/(2\sqrt b)$.
Each identity assumes the derivative and the quantities in its formula are defined.

| Rule | Expression | Derivative |
| --- | --- | --- |
| D05 | `Absolute(b)` | $\operatorname{Re}(\bar b \, \delta b) / \lvert b \rvert$, with `real_part`, so the derivative is real |
| D06 | `SquareRoot(b)` | $\delta b / (2 \sqrt b)$ |
| D07 | `Logarithm(b)` | $\delta b / b$ |

### Determinant and inverse

The determinant can have a nonzero derivative at a singular matrix. For
$A(t)=\operatorname{diag}(t,1)$, the determinant is $t$, so its derivative
at $t=0$ is one. Cofactors give this derivative directly. The inverse
rules have their own rank requirements.

| Rule | Expression | Derivative |
| --- | --- | --- |
| D08 | `Determinant(A)`, shape `(n, n)` | $\sum_{i,j=0}^{n-1} C_{ij}(A)\,\delta A_{ij}$, including at singular matrices |
| D09 | `Inverse(A)`, equal axis extents, and $A$ invertible | `-inverse(A) @ delta_A @ inverse(A)`, with `delta_A` = $\delta A$ |
| D10 | `Inverse(A)`, $A$ real with shape `(n, k)`, $k<n$, and full column rank | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |

In D08, $C_{ij}(A)=(-1)^{i+j}\det(A_{\widehat i,\widehat j})$ is the cofactor
obtained by deleting row $i$ and column $j$. The derivative must be valid for
every square matrix, including singular ones. For `n >= 2`, a rewrite may use
the cofactor sum or, by multilinearity, the sum of `n` determinants. In term
$i$ of the latter sum, replace row $i$ of $A$ with row $i$ of $\delta A$ and
keep the other rows of $A$. For a directional derivative, evaluate those
other rows at the differentiation point. For `n == 1`, the result is
`delta_A[0, 0]`. For `n == 0`, it is zero with the type of the derivative.

## Nested derivative examples

Displayed results omit zero terms and use schematic fresh names;
rewriting need not simplify algebra.

Let `u`, `a`, `b`, `h` and `k` be real scalar functions in a constant space,
with `a`, `b`, `h` and `k` independent of `u`.

| Input | Rule trace and result |
| ----- | --------------------- |
| `directional_derivative(directional_derivative(u * u, u, h), u, a, k)` | The inner derivative gives `h * u + u * h` by D01, D11 and SB18. Its point `u` is free under SB04. The outer derivative gives `h * k + k * h`. |
| `directional_derivative(directional_derivative(u * u, u, b, h), u, a, k)` | The inner derivative gives `h * b + b * h`. SH15 then finds no free `u`, so L10 gives the typed zero. |
