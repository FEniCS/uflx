# Evaluation and lowering

## Evaluation

A rewrite replaces an expression by one with the same mathematical value,
[binding](#substitution-and-binding) and [compatible type](#type-compatibility-of-rewrites).
It need not preserve floating-point results or simplify algebra. The
[mathematical admissibility](rules.md#mathematical-admissibility) convention
applies: zero times an undefined value remains undefined.

The lazy classes are `Action`, `DirectionalDerivative`, mesh `Integral`,
`PullBack`, `Gradient`, `Divergence`, `Curl` and `PushForward`.

The following operations commute through the listed operands when all other
operands are independent of the operation's target. Specific rules take
precedence. Make binders fresh for inserted operands before commuting.

| Node class | Operands | Operations that commute |
| ---------- | -------- | ----------------------- |
| `Sum`, `Tensor` | Every operand | `Action`, `DirectionalDerivative`, spatial partial derivative |
| `Indexed`, `ComponentTensor`, `IndexSum`, `Form`, `DirectionalDerivative` | The expression | `Action`, `DirectionalDerivative`, spatial partial derivative |
| `Division` | The numerator | `Action`, `DirectionalDerivative`, spatial partial derivative |
| `Gradient`, `Divergence`, `Curl`, `PushForward`, `PullBack`, `Integral` | The expression | `Action`, `DirectionalDerivative` |
| `Conjugate`, `RealPart`, `ImaginaryPart` | The expression | `DirectionalDerivative`, spatial partial derivative |

Directional differentiation uses a real parameter and holds domains fixed.
`Action` has its own [conjugation rule](#applying-actions); `RealPart` and
`ImaginaryPart` have no argument spaces. Actions on gradients and pull-backs
stay outside them, as specified in their rule tables.

### Transformations, recipes and contracts

| Component | Responsibility |
| --------- | -------------- |
| Rule | Local match, conditions and directed replacement; shared binding, type completion and reconstruction |
| Strategy | Rule selection and precedence, traversal, revisiting and stopping |
| Recipe | Composition of strategies with input, output and preservation contracts; the named passes below are reference recipes |

Rules do not schedule traversal or invoke a whole pass for cleanup. No
rule-object API or global registry is prescribed.

| Strategy condition | Requirement |
| ------------------ | ----------- |
| Specific identity or retained boundary | Takes precedence over generic expansion; omitting its rule does not enable an excluded fallback |
| Unmatched node | May remain during partial rewriting; required lowering must satisfy its output contract |
| Optional simplification | Need not eliminate any node class |
| Repetition | Requires progress; no cycles between opposite rewrites; a resource limit does not establish success |
| Structural stability | Does not imply mathematical canonicalisation or equal outputs from different recipes |

Contracts are independent predicates over the evaluation traversal:

| Contract | Required representation |
| -------- | ----------------------- |
| Core expressions | Structurally valid concrete core nodes and downstream terminals |
| Actions expanded | Every `Action` applies a `Function`, downstream terminal, `Gradient`, `PullBack` or another retained `Action` |
| Derivatives expanded | No `DirectionalDerivative` |
| Mesh integrals expanded | No `Integral` over a `Mesh` |
| Gradients expanded | Every `Gradient`, `Divergence` and `Curl` satisfies its retained-pattern condition in [The core normal form](#the-core-normal-form) |
| Pull-backs expanded | Every `PullBack` satisfies its retained-pattern condition in [The core normal form](#the-core-normal-form) |
| Reference integration | Derivatives expanded and mesh integrals expanded |
| Core normal form | The allowed patterns in [The core normal form](#the-core-normal-form) |

Contracts compose independently of schedules. Checks report offending nodes
and reasons without rewriting; they establish neither preservation of input
structure nor backend compatibility. Recipes separately state what domains or
operators they preserve. Quadrature and tabulation impose their own contracts.

### Standard lowering recipe

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `lower(f)` | | `f: AbstractExpression` | An equivalent expression in core normal form |

| `lower` condition | Requirement or result |
| ----------------- | --------------------- |
| Input | Core expressions and [prepared geometry](#evaluation-traversal); extension nodes need an extension recipe first |
| Unsupported rewrite | `TypeError` |
| Output violates core normal form | `ValueError`, identifying offending nodes |

The reference schedule is:

| Recipe step | Required result |
| ----------- | --------------- |
| [`apply_actions`](#applying-actions) | Actions expanded |
| [`apply_derivatives`](#applying-derivatives), including action cleanup | Actions and derivatives expanded |
| [`expand_integrals`](#expanding-integrals) | Also mesh integrals expanded |
| Repeat [`apply_gradients`](#applying-gradients), then [`expand_pull_backs`](#expanding-pull-backs), with action cleanup, until neither changes the expression | Core normal form |

This order evaluates gradients on their stated domains before transport:
full-dimensional `gradient(coordinate(D), D)` becomes the identity before
introducing Jacobians. Piola patterns survive until a matching cell pull-back.

Other recipes may select individual rules or families, subject to their local
conditions and geometry preparation contract. For example,
[shape differentiation](extensions.md#shape-derivatives) lowers spatially before
and after differentiation. Index-notation simplification, such as reducing an
`Indexed` of a `ComponentTensor`, remains a compiler operation.

## Rewriting and substitution

### Type compatibility of rewrites

A replacement $r$ must preserve the value of $e$ on its original domain and
satisfy these conditions:

| Attribute | Compatibility |
| --------- | ------------- |
| Value shape, scalar type, ordered argument spaces | Equal |
| Free indices | Same names and extents, up to order and consistent bound-name renaming |
| Domain | Equal, or $r$ has no domain and represents a constant field on $e$'s domain for each free-index and argument assignment |

For example, a constant pull-back or the gradient of coordinates may lose its
domain. A consumer needing that domain must retain it.

Before returning a table's result `r` or rebuilding its parent, complete its
type against the replaced expression `e`. Let

```python
z = Zero(e.value_shape, e.free_indices, e.argument_spaces, e.scalar_type)
```

Apply both steps in order, using tensor arithmetic:

| Condition | Completion |
| --------- | ---------- |
| `e` complex, `r` real | `Complex(Integer(1), Integer(0)) * r` |
| `r` lacks a free index of `e` | `r + z` |

The raw result must already have the required shape and argument spaces,
matching extents for retained indices and no extra free indices. A real `e`
requires a real result. A row giving zero uses `z`; this is what the tables
mean by zero with the type of an expression.

`Action` and `DirectionalDerivative` also promote real values before
substitution for complex functions, member by member for tuples; `replace`
does not. Promotion must precede evaluation of, for example, a square root at
a negative real point. Index completion must precede rebuilding a binder.

### Evaluation traversal

| Dependency | Evaluation |
| ---------- | ---------- |
| Root and value operands | Visit recursively in construction order; tuples member by member, including both `Action` operands |
| Binding declarations; non-expression `Base` objects | Do not visit |
| Domains, maps, spaces, indices and labels | May inspect, but do not traverse or modify stored expressions |
| Expressions exposed by a rule | Visit as separate roots; includes newly built map parameterisations |

Reference recipes require exposed geometry to satisfy core expressions,
actions expanded, derivatives expanded and mesh integrals expanded. The supplier
prepares it if necessary. This is a recipe precondition, not a constructor
restriction; another strategy may process unprepared geometry with suitable rules.

Normal-form claims cover only this traversal. Stored dependencies are checked
separately when needed; labels have no lowering requirement. Substitution uses
its broader dependency and binding rules below.

### Substitution and binding

Value operands are `AbstractExpression` construction arguments, including
those in tuples. Terminals have none. Labels are identity metadata, handled by
[Derived indices](expressions.md#derived-indices). Only these nodes bind names:

| Binding node | Value operands | Bound names | Bound in |
| ------------ | -------------- | ----------- | -------- |
| [`ComponentTensor(e, indices)`, `IndexSum(e, indices)`](expressions.md#componenttensor-and-indexsum) | `e` | The indices of `indices` | `e` and `indices` |
| [`Form(e, arguments)`](forms.md#form) | `e` | The functions of `arguments` | `e` and `arguments` |
| [`DirectionalDerivative(f, x, a, h)`](geometry.md#directionalderivative) | `f`, `a`, `h` | $x$ | $f$ and $x$ |

Tuple-valued operands are visited member by member. Substitution follows
construction dependencies subject to these scopes and the label rules;
individual passes specify which value operands they transform.

A binding node is made fresh for expressions $s_1, \dots, s_n$ by renaming,
in the operands in which the names are bound, each bound name that is free in
some $s_k$, so that the binding node does not capture the name.

Substitution replaces $r_1, \dots, r_n$ by $s_1, \dots, s_n$
simultaneously. Within a binding scope:

1. Drop each pair whose key $r_k$ has a bound name free in it.
2. Make the binding node fresh for the remaining replacement values $s_k$.
3. Substitute the remaining pairs into the operands in that scope.

Operands outside the binding scope use the original pairs. Inserted values
are never searched for further matches.

For freshening against $s_1, \ldots, s_n$, let `body` be the tuple of operands
in the binder's scope. Use these names:

| Bound name | Fresh name |
| ---------- | ---------- |
| Index `i` | `Index(("uflx.rename", i, body, (s_1, ..., s_n)))` |
| Function `v` in $V$ | `Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))` |

The [derived labels](expressions.md#derived-indices) keep these names distinct
and absent from `body` and the inserted expressions.

The core pass `replace(e, mapping)` performs this simultaneous substitution,
replacing each key of the `dict` `mapping` with its value.

Match whole original nodes first, subject to binding scope. A match stops
descent; inserted values and rebuilt nodes are not searched. Otherwise,
recurse according to the binding and label rules. Overlapping keys are
allowed: the outermost eligible match wins.

Order replacement pairs by the structural key order of
[`Base`](rules.md#the-base-object-base), preserving that order through scope filtering
and in fresh-name labels. Dictionary insertion order has no effect.

```python
replace(a + b, {a + b: c, a: d})  # c
replace(a + b, {a: d, d + b: c})  # d + b; no match after rebuilding
replace(a, {a: b, b: c})  # b; inserted values are not searched
```

Examples, with `j` a fresh index:

| Substitution | Result |
| ------------ | ------ |
| `replace(index_sum(a[i] * c, (i,)), {c: b[i]})` | `index_sum(a[j] * b[i], (j,))`; `i` remains free |
| `replace(form(a, (v,)), {v: w})` | The original form; `v` is bound |
| `replace(DirectionalDerivative(f, u, u, h), {u: u0})` | `DirectionalDerivative(f, u, u0, h)`; the point is outside the binding scope |

Substitution follows all construction dependencies, including domains, maps,
spaces and derived labels. It commutes with a fixed composition when the
static-type choice is unchanged: replacing `A` by `B` in `transpose(A)`
gives `transpose(B)`, with the derived indices of `B`.

### Grammars

A grammar is a set of expressions describing a structural contract. A check
asks whether an expression belongs to that set (see
[New questions](extensions.md#new-questions)). Contracts refer only to nodes
visited by the [evaluation traversal](#evaluation-traversal).

A grammar pattern inspects a node and objects within a bounded number of edges,
not the result of a pass. It may classify a domain or map without adding that
dependency to the traversal. For a grammar $N$ and patterns $P$, $N \setminus P$
contains the expressions of $N$ with no visited node matching $P$. Contracts
compose by intersection; there is no mandatory sequence of numbered grammars.

A recipe maps an input grammar to an output grammar under its stated
preconditions. A single rule only guarantees its replacement and preservation
obligations. Extension nodes require extension grammars and rules before a
consumer can require core expressions.

The form layer consists of every `Action` and every node with nonempty argument
spaces; the integrand layer consists of the other nodes. A contract may restrict
these separately. A `Function` can satisfy core normal form even if its stored
parameterisation contains a reducible gradient: that parameterisation is checked
separately when exposed.

## Actions and derivatives

### Applying actions

| `apply_actions(f)` |  |
| ------------------ | - |
| Evaluates          | `Action` |
| Input contract     | Core expressions |
| Output contract    | Actions expanded |
| Strategy | Inner actions first; first matching row; revisit created actions, including those exposed by substitution |

Downstream terminals may declare argument spaces, for example for assembled
vectors or matrices; their libraries evaluate the retained actions.

| Form $f$ | Rewrite of `Action(f, p, g)` |
| -------- | ---------------------------- |
| `Form(e, arguments)`, $p$ at least the arity $n$ of $e$ | The value of `Action`: `Form(e, others)` with $g$ substituted for the argument at $p - n$, or $e$ with $g$ substituted if no entries remain |
| `Conjugate(e)` | `conjugate(action(e, p, conjugate(g)))`, with conjugation member by member for a `tuple` $g$ |
| `Product(a, b)`, with $n$ the arity of $a$ | `Product(Action(a, p, g), b)` for $p < n$, otherwise `Product(a, Action(b, p - n, g))` |
| `Gradient`, `PullBack` | Remains; the operation transforms the form's output |
| A commuting node class (see [Evaluation](#evaluation)) | The node with the action moved into the listed operands |
| `Zero` | Zero with the type of the action |
| `Function` or a downstream terminal | Remains, for the supplying library |
| `Action` | Remains, as a further application of a remaining action |
| Any other class of the core | Cannot occur, because the argument spaces are `()` |

Substitution preserves all enclosing and inserted binders. In particular,
commuting through `DirectionalDerivative` makes its variable fresh for `g`;
occurrences in `g` are not differentiated.

### Derivatives of the node classes

Both differentiation contexts use this table after their specific rules.
$\delta b$ is the operand's derivative; each context defines independence from
its differentiation target. Extensions add derivative rules here for reuse.

| Expression | Derivative |
| ---------- | ---------- |
| `Product(b, c)` | $\delta b \, c + b \, \delta c$ |
| `Division(b, c)` | $\delta b / c - (b \, \delta c) / (c c)$ |
| `Power(b, c)`, $c$ constant | $c \, b^{c - 1} \, \delta b$ |
| `Power(b, c)` | $b^c (\log(b) \, \delta c + c \, \delta b / b)$ |
| `Absolute(b)` | $\operatorname{Re}(\bar b \, \delta b) / \lvert b \rvert$, with `real_part`, so the derivative is real |
| `SquareRoot(b)` | $\delta b / (2 \sqrt b)$ |
| `Logarithm(b)` | $\delta b / b$ |
| `Determinant(A)`, $A$ invertible | $\det(A)\sum_{i,j}(A^{-1})_{ij}\,\delta A_{ji}$, with `inverse`, `index_sum` and derived indices $i,j$ |
| `Inverse(A)`, equal axis extents, including symbolic extents | `-inverse(A) @ delta_A @ inverse(A)`, with `delta_A` = $\delta A$ |
| `Inverse(A)`, rectangular shape | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |
| `Conjugate(b)`, `RealPart(b)`, `ImaginaryPart(b)` | The same node applied to $\delta b$ |
| A commuting node class (see [Evaluation](#evaluation)) | The node with the derivative moved into the listed operands |

The constant-exponent power rule avoids a logarithm, permitting negative real
bases where differentiable. Results are core expressions; recipes evaluate
resulting actions using [Applying actions](#applying-actions).

### Applying derivatives

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` and actions created by its evaluation |
| Input contract         | Core expressions; actions expanded |
| Output contract        | Actions and derivatives expanded |
| Strategy | Inner derivatives first; first matching row; alternate with `apply_actions` until both output contracts hold |

Write $\mathrm{D}b$ for `DirectionalDerivative(b, x, a, h)`. Match original
operands before substituting `[x := a]` into undifferentiated operands; insert
`h` unchanged. Rules may produce actions; cleanup belongs to the recipe.

| Expression $f$ | Rewrite of the derivative at $a$ |
| -------------- | -------------------------------- |
| No function of $x$ free in $f$ | Zero with the type of the derivative |
| $f$ equal to the member $x_k$ of $x$ | $h_k$, or $h$ for a function $x$ |
| `Action(b, p, g)` | `Action(Db, p, g) + Action(b, p, Dg)`, with undifferentiated operands evaluated at $a$ |
| `DirectionalDerivative` | Defer until the inner derivative is expanded |
| Any other expression | The row of [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \mathrm{D} b$ |

The dependency test follows value operands and binding scopes, holding domains
and maps fixed. A `Function` contributes itself, not names in its label or
space; a form binding `x` has no free `x`. Tuple arguments vary member by member,
including conjugated entries. Dual directions and substitution points may be
symbolic forms whose actions the recipe evaluates.

For real scalar functions `u` and `a` in a space of constants,
`directional_derivative(u * u, u, a, u)` rewrites to `u * a + a * u`.

For a real scalar function `u` in a space of constants and a constant array
`a` of shape `(m,)`, `directional_derivative(u, u, a[i], Integer(1))` lowers
to `Integer(1) + Zero((), ((i, m),), (), real)`: its value is one for every `i`.
Type completion retains indices even from unused members of tuple points
and directions.

## Spatial lowering

### Expanding integrals

| `expand_integrals(f)` |  |
| --------------------- | - |
| Evaluates             | `Integral` over a mesh |
| Input contract        | Core expressions |
| Output contract       | Mesh integrals expanded |
| Strategy | Inner mesh integrals first |

For group $j$ of `D.topology`, let $K_j$ be the reference cell and $c_j$ =
`cell_index(D.topology, j)`.

Choose the internal index without capturing $f$: form the template
`index_sum(pull_back(Integer(1), cell_map(D, j)), (c_j,))`, make its binder
fresh for $f$ using [Substitution and binding](#substitution-and-binding), and
let $e_j$ be its resulting bound index. This retains $c_j$ when no renaming is
needed, preserving shared geometry for ordinary integrands. Set
`D_j = cell_map(D, j, e_j)` and
`piece_j = domain_map(K_j, D_j.parameterisation, D.ambient_domain)`.

Let `K` contain one derived `Index` for each value axis of $f$, with `()`
for a scalar, and set `f_K = indexed(f, K)`.

| Integral | Rewrite |
| -------- | ------- |
| `Integral(f, D)` | `component_tensor(sum_j, K)`, where `sum_j` is the sum over groups $j$ of `index_sum(integral(g_j * density(D_j), K_j), (e_j,))` |

$g_j$ is `pull_back(f_K, piece_j)` when $f$ has domain `D.ambient_domain`,
and `pull_back(f_K, D_j)` otherwise, including domain-free $f$.
The pull-back introduces the internal cell index; `IndexSum` consumes it and
retains every free index of $f$ and the component indices `K`. The outer
`component_tensor` restores the value shape, and is the identity for scalar
$f$. The remaining integrals have coordinate domains.

Repeated expansion of the same input constructs the same graph. Paths involving
capture-avoiding renaming may choose different bound names; the rule does not
assert structural equality between every ordering of substitution and expansion.
They must preserve the same value and free indices.

### Applying gradients

| `apply_gradients(f)` |  |
| -------------------- | - |
| Evaluates            | `Gradient`, `Divergence` and `Curl` |
| Input contract       | Core expressions; actions and derivatives expanded |
| Output contract      | Actions, derivatives and gradients expanded |
| Strategy | Inner expressions first; first matching row; evaluate exposed actions with `apply_actions` and revisit affected expressions |

A cell is a `CellMap` of a topology or mesh. The divergence and curl boundary
rows identify the Piola patterns. Families select rules, not public passes.

$\partial_j b$ denotes component $j$ of the gradient of $b$ on $P$, with
$b$'s value shape and the additional free index $j$ of ambient extent $n$.
Gradient reconstruction uses derived component indices $A$ and $j$, one
index in $A$ per value axis of $f$. Empty lists use the factory identities.

| Family | Expression or operation | Rewrite |
| ------ | ----------------------- | ------- |
| Action transport | `Gradient(Action(B, p, g), P)` | `Action(Gradient(B, P), p, g)` |
| Boundary | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)` | Remains for the Piola rule of a cell pull-back |
| Boundary | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$ | Remains for the Piola rule of a cell pull-back |
| Differential expansion | Other `Divergence` or `Curl` | Its composition of gradients from [Divergence and Curl](geometry.md#divergence-and-curl) |
| Gradient reconstruction | Other `Gradient(f, P)` | `component_tensor(∂_j f[A], A + (j,))`, using the component rules below |
| Gradient rules | $\partial_j b$, $b$ with no domain | Zero with $b$'s shape and argument spaces, and its free indices followed by $j$ |
| Gradient rules | $\partial_j$ `Coordinate(P)`, $P$ of its ambient dimension | `component_tensor(I[i, j], (i,))`, with $I$ the identity of shape `(n, n)` and a derived $i$ |
| Boundary | $\partial_j b$, $b$ a terminal, `Gradient`, `PushForward`, Piola-matched node, or `PullBack` by a cell or `SideMap` | `component_tensor(gradient(b, P)[B + (j,)], B)`, with one derived index in `B` per value axis of $b$; the gradient remains |
| Gradient rules | $\partial_j$ `PullBack(g, φ)`, φ a `DomainMap` | `component_tensor(index_sum(H[A + (l,)] * J[l, j], (l,)), A)`, with `H = pull_back(gradient(g, Q), φ)`, `J = jacobian(φ)`, and $Q$ the ambient domain of φ |
| Action transport | $\partial_j$ `Action(B, p, g)` | `Action(∂_j B, p, g)`, with $g$ unchanged |
| Gradient rules | $\partial_j b$, any other $b$ | [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \partial_j b$ |

In the chain rule, `A` indexes the value axes of `g` and `l` is a distinct
contraction index; scalar `g` uses `A == ()`. On embedded meshes the same rules
apply to the tangential gradient: the coordinate gradient remains a projection
until cell transport. Retained terminal gradients belong to the supplying
library; gradients exposed by pull-backs are processed in the next spatial round.

### Expanding pull-backs

| `expand_pull_backs(f)` |  |
| ---------------------- | - |
| Evaluates              | `PullBack` |
| Input contract         | Core expressions; actions and derivatives expanded; prepared geometry when exposed |
| Output contract        | Actions, derivatives and pull-backs expanded |
| Strategy | Inner pull-backs first; first matching row; revisit created pull-backs and evaluate exposed actions with `apply_actions` |

$P$ and $Q$ are the parameter domain and the ambient domain
of $\psi$. For $\psi$ a cell $D_j$ = `cell_map(D, j, e)` of a mesh $D$ with
the topology $T$, $\hat f$ is `pull_back(f, cell_map(T, j, e))` and $m_\psi$
is the law $m$ on $\psi$ (see [`TransformationLaw`](geometry.md#transformationlaw)). For
$\psi$ = `cell_map(Y, i, q)` for a topology or mesh $Y$ and a side $s$,
`B_selected` and
`A_selected` are the maps of [Restriction](geometry.md#restriction), selected by `q`.

| Family | `PullBack(f, ψ)` | Rewrite |
| ------ | ---------------- | ------- |
| Action transport | `Action(A, p, g)` | `Action(PullBack(A, ψ), p, g)`, with $g$ unchanged |
| Binding | `Form(e, arguments)` | Make the form fresh for $\psi$, obtaining `Form(e_new, arguments_new)`; return `Form(PullBack(e_new, ψ), arguments_new)`, preserving argument spaces |
| Boundary | $\psi$ a `SideMap`, $f$ with a domain | Remains until a pull-back by a cell of $P$ |
| Pointwise | $f$ with no domain | $f$, with [type completion](#type-compatibility-of-rewrites) to retain the free indices of $\psi$ |
| Geometry | `Coordinate(Q)` | The parameterisation of $\psi$ |
| Restriction | `PullBack(g, SideMap(X, Y, s))`, $\psi$ = `cell_map(Y, i, q)` | `PullBack(PullBack(g, B_selected), A_selected)` |
| Boundary | `PullBack(g, φ)`, φ a `SideMap` or a cell | Remains, so a pull-back nests only inside a read on a cell |
| Map composition | `PullBack(g, φ)`, φ a `DomainMap` | `PullBack(g, domain_map(P, pull_back(φ.parameterisation, ψ), φ.ambient_domain))`, the pull-back by the image of $\psi$ under φ (see [`PullBack`](geometry.md#pullback)) |
| Mapped values | `PushForward(f, D, m)`, $\psi$ a cell of $D$ | $m_\psi(\hat f)$ |
| Piola | `Divergence(PushForward(f, D, m), D)` of the identity of Piola, $\psi$ a cell of $D$ | $m'_\psi$ of `divergence(f̂, K_j)`, with $m'$ = `transformation_law((), 1)` |
| Piola | `Curl(PushForward(f, D, m), D)` of the identity of Piola, $\psi$ a cell of $D$ | $m'_\psi$ of `curl(f̂, K_j)`, with $m'$ = `transformation_law(("contravariant",), 1)` for $n = 3$ and `transformation_law((), 1)` for $n = 2$ |
| Boundary | `PushForward`, `Divergence` or `Curl` | Remains, read at the points of $P$ |
| Gradient transport | `Gradient(g, Q)`, $\psi$ a cell of $Q$ | $m_\psi(\nabla_P(\psi^* g))$, with $m$ = `transformation_law(("covariant",), 0)` |
| Boundary | `Gradient(g, Q)` | Remains, the gradient on $Q$ read at the points of $P$ |
| Geometry | `Normal(SideMap(D, E, s))`, $\psi$ = `cell_map(E, i, q)` | $w / \lvert w \rvert$ of [`Normal`](geometry.md#normal) |
| Boundary | A terminal | Remains, for the library of the terminal |
| Pointwise | Any other node class | Make its binders fresh for $\psi$; rebuild from the pull-backs of its expression operands |

Map composition processes its new parameterisation as a separate root under
the [geometry preparation contract](#evaluation-traversal). Cell maps retain
explicit selection; unselected side restrictions remain. Piola rows read no
derivative of $J$. Only cell transport inverts a Jacobian; other maps use the
[gradient chain rule](#applying-gradients). Libraries evaluate remaining cell
reads under [Functions in a space](forms.md#functions-in-a-space).

## The core normal form

Core normal form is the target of the reference lowering recipe. It combines
core expressions, actions expanded, derivatives expanded, mesh integrals
expanded, gradients expanded and pull-backs expanded. Its expression nodes
are visited by the evaluation traversal.
The following table lists all allowed node classes and their conditions:

| Chapter | Node class | Allowed in the core normal form |
| ------- | ---------- | ------------------------------- |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Integer`, `Float`, `Complex`, `Zero` | Always |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Sum`, `Product`, `Division`, `Power`, `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm` | Always |
| [Tensors and indexing](expressions.md#tensors-and-indexing) | `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | Always, including index notation that a rewrite could reduce, such as an `Indexed` of a `ComponentTensor`, which the compiler rewrites (see [Evaluation](#evaluation)) |
| [Tensor operations](expressions.md#tensor-operations) | `Inverse`, `Determinant` | Always, for the compiler to evaluate |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Coordinate` | Always |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Gradient` | Of a terminal, a `Gradient`, a `PushForward`, a node of the identity of Piola, or a `PullBack` by a cell or a `SideMap`, other than `Coordinate(P)` on its own domain $P$ of the dimension of its ambient coordinate domain |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Divergence`, `Curl` | Of the identity of Piola, under no pull-back by a cell |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Integral` | Over a domain that is not a `Mesh`, such as a reference cell or `RD(0)` |
| [Maps between domains](geometry.md#maps-between-domains) | `PullBack` | By a `SideMap` of an expression other than an `Action` or `Form`, or of a terminal, of a `PullBack` by a `SideMap` or a cell, or of a `Gradient`, a `PushForward` or a node of the identity of Piola by a map that is not a cell |
| [Meshes](geometry.md#meshes) | `Normal` | Always |
| [Meshes](geometry.md#meshes) | `PushForward` | Under no pull-back by a cell |
| [Function spaces](forms.md#function-spaces) | `Function` | Always |
| [Variational forms](forms.md#variational-forms) | `Form` | Always, including a form whose expression is a form (see [`Form`](forms.md#form)) |
| [Variational forms](forms.md#variational-forms) | `Action` | Its form is a `Function`, a downstream terminal, a `Gradient`, a `PullBack` or another such `Action`, with each operand in the core normal form |
| | A terminal of a downstream library | Always |

Dependencies outside the [evaluation traversal](#evaluation-traversal), such
as geometry and labels, are outside this normal-form claim.

The code that assembles receives the form layer (see [Grammars](#grammars)):
sums and products of forms, real and complex zeros, forms whose expression
may be a form, with the order of the arguments from the static type (see [`Form`](forms.md#form)),
and actions on supplied or spatially transformed forms.

A form compiler receives the integrand layer: integrals over reference cells
and points, with arithmetic and index notation over constants, reference
coordinates and the library reads of [Functions in a space](forms.md#functions-in-a-space).

```{dropdown} Reusing geometry calculations
Every expression that a pass builds from domain maps, such as the Jacobian,
the inverse Jacobian, the density and the normal of a side, depends only on
those domain maps and remains a subexpression of the result. A form compiler
can therefore evaluate such an expression once for the values that the domain
maps read, such as the coordinates of a mesh, and use the value in every
integral that contains the expression.
```
