# Applications and mappings

These reference recipes use the [shared rewrite procedure](rewriting.md#visiting-and-revisiting)
and [representation contracts](evaluation.md#representation-contracts).

## Expanding applications

| `expand_applications(f)` |  |
| ------------------ | - |
| Evaluates          | `Apply` |
| Input contract     | Core expressions |
| Output contract    | Applications expanded |

| Selection | Rules |
| --------- | ----- |
| Apply explicit forms | SB16 |
| Conjugation and products | L05, ST02 |
| Keep spatially transformed forms | ST03 |
| Commute into operands | L01, ST01, L02–L03 |
| Zero and remaining applications | L07, ST04, ST06 |

### Application rules

An application supplies every argument of a form. Substitution and commuting
preserve binding under SB01–SB15 and SB20–SB21. The separate `partial` operation
uses the retained-name rejection rule in [Partial application](forms.md#partial-application).

| Rule | Form $f$ | Rewrite of `Apply(f, G)` |
| --- | --- | --- |
| SB16 | `Form(e, entries)`, with $n=$ `len(e.argument_spaces)` | In `Apply(Form(e, entries), G)`, substitute `G[n:]` for `entries` simultaneously in `e`, yielding `e_new`, then return `apply(e_new, G[:n])`. For a conjugated entry substitute the conjugate value for its function; handle product-space members simultaneously. Promote before substitution and preserve binding by SB07–SB14. Inserted values are not substituted into one another. |
| ST02 | `Product(a, b)`, with $n$ the arity of $a$ | `Product(apply(a, G[:n]), apply(b, G[n:]))` |
| ST03 | `Gradient`, `PullBack` | Remains; the operation transforms the form's output |
| ST04 | `Function` or a downstream terminal | Remains, for the supplying library |
| ST06 | An expression with `argument_spaces == ()` | `apply(f, ())` returns `f`; an `Apply` node requires at least one slot |
| L05 | `Conjugate(e)` | `conjugate(apply(e, conjugated_G))`, where `conjugated_G` conjugates every supplied member, preserving the tuple grouping of `G` |
| L07 | `Zero` | Zero with the type of the application |

| Rule | Applies to | Requirement |
| --- | --- | --- |
| SB17 | Binding during application | Preserve enclosing and inserted binders. Before commuting an application through a derivative or index binder, make its bound names fresh for every supplied value; inserted values are not captured or differentiated by that inner operation. |

For commuting through other nodes, use the shared [structural](rewriting.md#structural-rules)
and [linearity](rewriting.md#linearity) rules. Downstream terminals may declare argument
spaces, for example for assembled vectors or matrices; their libraries evaluate
remaining applications.

## Expanding integrals

| `expand_integrals(f)` |  |
| --------------------- | - |
| Evaluates             | `Integral` over a mesh |
| Input contract        | Core expressions |
| Output contract       | Mesh integrals expanded |

| Selection | Rules |
| --------- | ----- |
| Mesh integration | I01 with the notation and index construction in [Integration](#integration) |

### Integration

Mesh integration produces a local contribution by adding the integrals over
its statically known cell slots. For slot $j$ of `D.topology`, set
`K_j = D.topology.cell_types[j]`, `D_j = cell_map(D, j)` and
`piece_j = domain_map(K_j, D_j.parameterisation, D.ambient_domain)`.

Let `K` contain one derived `Index` for each value axis of `f`, with `()`
for a scalar, and set `f_K = indexed(f, K)`.

| Rule | Integral | Rewrite |
| ---- | -------- | ------- |
| I01 | `Integral(f, D)`, $D$ a `Mesh` | `component_tensor(sum_j, K)`, where `sum_j` adds `integral(g_j * density(D_j), K_j)` over the fixed local slots $j$, in tuple order |

$g_j$ is `pull_back(f_K, piece_j)` when $f$ has domain `D.ambient_domain`,
and `pull_back(f_K, D_j)` otherwise, including domain-free $f$.
No cell index is created or bound. The outer `component_tensor` restores the
value shape and is the identity for a scalar. All free indices of `f` remain
free. The remaining integrals have coordinate domains. The package performs
mesh iteration and assembly outside this rewrite.

## Expanding pull-backs

| `expand_pull_backs(f)` |  |
| ---------------------- | - |
| Evaluates              | `PullBack` |
| Input contract         | Core expressions; applications and derivatives expanded; prepared geometry when exposed |
| Output contract        | Applications, derivatives and pull-backs expanded |

Run `expand_applications` on exposed applications and continue under the common
procedure. Leave new gradients to a later call to `apply_gradients`.

| Selection | Rules |
| --------- | ----- |
| Applications and binding | ST11, SB19 |
| Keep a side restriction | R01 |
| Constants and coordinates | PB01–PB02 |
| Select the cell on a side | R02 |
| Nested maps | PB03–PB04 |
| Mapped values | CC02 |
| Piola identities | CC08–CC09 |
| Remaining mapped derivatives | CC03–CC05 |
| Side normal | R03 |
| Terminals and pointwise operations | PB05–PB06 |

### Applications and forms under pull-back

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| ST11 | `Apply(A, G)` | `Apply(PullBack(A, ψ), G)`, with every supplied value unchanged |
| SB19 | `Form(e, arguments)` | Make the form fresh for $\psi$, obtaining `Form(e_new, arguments_new)`; return `Form(PullBack(e_new, ψ), arguments_new)`, preserving argument spaces |

### Pull-back algebra

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| PB01 | $f$ with no domain, neither an `Apply` nor a `Form` | The same spatially constant value $f$ under [Domain-free expressions](expressions.md#domain-free-expressions), with [type completion](rewriting.md#type-compatibility-of-rewrites) to retain the free indices of $\psi$ |
| PB02 | `Coordinate(Q)`; exclude R01 | The parameterisation of $\psi$ |
| PB03 | `PullBack(g, φ)`, φ a `SideMap` or a cell; exclude R01 and R02 | Remains, so a pull-back nests only inside a read on a cell |
| PB04 | `PullBack(g, φ)`, φ a `DomainMap`; exclude R01 | `PullBack(g, domain_map(P, pull_back(φ.parameterisation, ψ), φ.ambient_domain))`, the pull-back by the image of $\psi$ under φ (see [`PullBack`](geometry.md#pullback)) |
| PB05 | A terminal not matching PB01, PB02, R01 or R03 | Remains, for the library of the terminal |
| PB06 | $f$ has a domain, $\psi$ is not a `SideMap`, and $f$ belongs to a class in the table below | Make its binders fresh for $\psi$; rebuild from the pull-backs of the listed operands, preserving the other construction arguments |

| PB06 operand class | Operands to pull back |
| ------------------ | --------------------- |
| `Sum`, `Product`, `Division`, `Power` | Both expression operands |
| `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm`, `Inverse`, `Determinant` | The expression operand |
| `Tensor` | Each entry |
| `BroadcastTo` | `value`; keep the target shape |
| `Indexed` | `expression`; keep the index entries |
| `ComponentTensor`, `IndexSum` | `expression`, after making the bound indices fresh for $\psi$ |

Other classes need their own rules; adding a node class does not extend PB06.

Map composition processes its new parameterisation as a separate root under
the [geometry preparation contract](rewriting.md#evaluation-traversal).

### Change of coordinates

For a mesh cell $\psi = D_j =$ `cell_map(D, j)` with topology $T$, let
$K_j$ be its reference cell, $\hat f$ = `pull_back(f, cell_map(T, j))`,
and $m_\psi$ be the law $m$ on $\psi$ (see
[`TransformationLaw`](geometry.md#transformationlaw)).

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| CC02 | `PushForward(f, D, m)`, $\psi$ a cell of $D$ | $m_\psi(\hat f)$ |
| CC03 | `PushForward`, `Divergence` or `Curl` not matching CC02, CC08 or CC09; exclude R01 | Remains, read at the points of $P$ |
| CC04 | `Gradient(g, Q)`, $\psi$ a cell of $Q$ | $m_\psi(\nabla_P(\psi^* g))$, with $m$ = `transformation_law(("covariant",), 0)` |
| CC05 | `Gradient(g, Q)` not matching CC04; exclude R01 | Remains, the gradient on $Q$ read at the points of $P$ |

Only cell transport inverts a Jacobian; other maps use the chain rule.
Libraries evaluate remaining cell reads under
[Functions in a space](forms.md#functions-in-a-space).

### Piola identities

These cell pull-back rules consume the [retained Piola patterns](differentiation.md#retained-piola-patterns)
without differentiating the Jacobian.

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| CC08 | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)`, $\psi$ a cell of $D$ | $m'_\psi$ of `divergence(f̂, K_j)`, with $m'$ = `transformation_law((), 1)` |
| CC09 | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$, $\psi$ a cell of $D$ | $m'_\psi$ of `curl(f̂, K_j)`, with $m'$ = `transformation_law(("contravariant",), 1)` for $n = 3$ and `transformation_law((), 1)` for $n = 2$ |

### Restriction

For $\psi =$ `cell_map(Y, i)` of a topology or mesh $Y$ and a side $s$,
`B_i` and `A_i` are the maps of
[Restriction](topology.md#restriction).
Cell maps retain explicit selection; unselected side restrictions remain.

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| R01 | $\psi$ a `SideMap`, $f$ with a domain and neither an `Apply` nor a `Form` | Remains until a pull-back by a cell of $P$ |
| R02 | `PullBack(g, SideMap(X, Y, entries))`, $\psi$ = `cell_map(Y, i)` | `PullBack(PullBack(g, B_i), A_i)` |
| R03 | `Normal(SideMap(D, E, entries), reference_conormals)`, $\psi$ = `cell_map(E, i)`; the cell Jacobian has full column rank and $w \ne 0$ | $w / \lvert w \rvert$, with $w$ defined in [The pull-back of Normal](geometry.md#the-pull-back-of-normal) |

### Retained pull-backs

For the pull-back representation contract, `PullBack` is allowed in exactly
these cases:

- By a `SideMap` of an expression other than an `Apply` or `Form`.
- Of a terminal.
- Of a `PullBack` by a `SideMap` or a cell.
- Of a `Gradient`, a `PushForward` or a node of a Piola identity, by a map
  that is not a cell.

These checks do not override rule exclusions or precedence.

## Expanding tensor arithmetic

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `expand_arithmetic(f)` | | `f: AbstractExpression` | An equivalent expression with tensor arithmetic expanded |

| `expand_arithmetic(f)` | |
| ---------------------- | - |
| Evaluates | Tensor-valued `Sum`, `Product` and `Division` |
| Input contract | Core expressions |
| Output contract | Tensor arithmetic expanded |

Apply ST13–ST16 to tensor arithmetic; retain scalar arithmetic.

### Tensor arithmetic components

For non-scalar results, `K` has one derived index per value axis, using
`Index(("uflx.expand_arithmetic", f, m))` for node `f` and axis `m`.
Choose fresh names under SB07–SB14; operand free indices remain free.

| Rule | Node | Rewrite |
| ---- | ---- | ------- |
| ST13 | `Sum(a, b)` | `component_tensor(Sum(a[K], b[K]), K)` |
| ST14 | `Product(a, b)`, scalar `a` | `component_tensor(Product(a, b[K]), K)` |
| ST15 | `Product(a, b)`, scalar `b` | `component_tensor(Product(a[K], b), K)` |
| ST16 | `Division(a, b)` | `component_tensor(Division(a[K], b), K)` |

## Piola example

Displayed results omit zero terms and use schematic fresh names;
rewriting need not simplify algebra.

### Piola identity exposed by an application

Let `v` and `u` be vector functions on the topology of mesh `D`, with geometry
independent of them. Set `m = transformation_law(("contravariant",), 1)`,
`F = form(push_forward(v, D, m), (v,))`, and let `ψ` be a cell map of `D`.

| Step | Result |
| ---- | ------ |
| Start with `pull_back(divergence(apply(F, (u,)), D), ψ)` | The Piola pattern is hidden by the application |
| Apply SB16 | `pull_back(divergence(push_forward(u, D, m), D), ψ)` |
| Apply CC06 | Retain the divergence until the cell pull-back is processed |
| Apply CC08 | The reference-cell divergence of `û`, carried by `transformation_law((), 1)` |

For reference field $\hat u=(\xi,\eta)$, the reference divergence is 2.

| Cell map | Physical divergence |
| -------- | ------------------- |
| $\psi(\xi,\eta)=(2\xi,3\eta)$ | $2/6=1/3$ |
| $\psi(\xi,\eta)=(\xi,(1+\xi)\eta)$, $\xi>-1$ | $2/(1+\xi)$ |

For the second map the physical field is
$u(x,y)=(x/(1+x),\;y(1+2x)/(1+x)^2)$. Direct differentiation gives
$1/(1+x)^2+(1+2x)/(1+x)^2=2/(1+x)$, agreeing with CC08 without derivatives
of the Jacobian in the lowered expression.

This order evaluates gradients on their stated domains before changing coordinates.
Full-dimensional `gradient(coordinate(D), D)` becomes the identity before
introducing Jacobians. Expressions matching Piola identities remain until a
matching cell pull-back.
