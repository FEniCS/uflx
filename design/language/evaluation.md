# Rewriting and lowering

Rewriting replaces expressions without changing their mathematical value.
Lowering applies rewrites to meet the requirements of later computations.

Rule IDs stay with their rules when moved. They cover rewrite rows, including
retained cases, and shared type, traversal and substitution rules. Conditions
and notation surrounding a table apply to its numbered rows. Contracts and
recipe schedules are not numbered. New rules take new IDs; existing IDs are
not renumbered or reused.

## Rewriting

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R001 | Every rewrite | A rewrite preserves [binding](#substitution-and-binding) and must have a [compatible type](#type-compatibility-of-rewrites). It need not preserve floating-point results or simplify algebra. The [mathematical admissibility](rules.md#mathematical-admissibility) convention applies. Zero times an undefined value remains undefined. |

The lazy nodes are `Action`, `DirectionalDerivative`, mesh `Integral`,
`PullBack`, `Gradient`, `Divergence`, `Curl` and `PushForward`. They store
operations to evaluate later.

### Transformations, recipes and contracts

| Component | Responsibility |
| --------- | -------------- |
| Rule | Match part of an expression, check conditions and replace it; use the shared rules for binding, type completion and rebuilding nodes |
| Strategy | Choose rules and their order, which nodes to visit, when to revisit them and when to stop |
| Recipe | Combine strategies with input, output and preservation requirements |

The named passes below are reference recipes. Rules cannot call a whole pass
to simplify their results. The spec requires no rule objects or global registry.

A contract states a requirement, such as “no directional derivatives remain”.
Contracts can be combined and apply to the nodes visited by the
[evaluation traversal](#evaluation-traversal).

| Strategy condition | Requirement |
| ------------------ | ----------- |
| Specific identity or a pattern that must remain | Takes precedence over generic expansion, even when its rule is not selected |
| Unmatched node | May remain during partial rewriting; required lowering must satisfy its output contract |
| Optional simplification | Need not eliminate any node class |
| Repetition | Must make progress without cycling; reaching a time or step limit does not count as success |
| Expression no longer changes | Does not mean all equivalent expressions have the same structure, or that different recipes give the same output |

## Rewriting and substitution

### Type compatibility of rewrites

A replacement $r$ must preserve the value of $e$ on its original domain and
satisfy these conditions.

| Rule | Attribute | Compatibility |
| ---- | --------- | ------------- |
| R002 | Value shape, scalar type, ordered argument spaces | Equal |
| R003 | Free indices | Same names and extents, up to order and consistent bound-name renaming |
| R004 | Domain | Equal, or $r$ has no domain and represents a constant field on $e$'s domain for each free-index and argument assignment |

For example, a constant pull-back or the gradient of coordinates may lose its
domain. Code that needs the original domain must keep it separately.

*Type completion* restores scalar type and free indices lost in a rewrite.
Apply it before returning a table's result `r` or rebuilding its parent.
For the replaced expression `e`, let

```python
z = Zero(e.value_shape, e.free_indices, e.argument_spaces, e.scalar_type)
```

Apply both steps in order, using tensor arithmetic.

| Rule | Condition | Completion |
| ---- | --------- | ---------- |
| R005 | `e` complex, `r` real | `Complex(Integer(1), Integer(0)) * r` |
| R006 | `r` lacks a free index of `e` | `r + z` |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R007 | Result before type completion | The raw result must already have the required shape and argument spaces, matching extents for retained indices and no extra free indices. A real `e` requires a real result. A row giving zero uses `z`; this is what the tables mean by zero with the type of an expression. |
| R008 | Promotion and completion order | `Action` and `DirectionalDerivative` also promote real values before substitution for complex functions, member by member for tuples; `replace` does not. Promotion must precede evaluation of, for example, a square root at a negative real point. Index completion must precede rebuilding a binder. |

### Evaluation traversal

| Rule | Dependency | Evaluation |
| ---- | ---------- | ---------- |
| R009 | Root and value operands | Visit recursively in construction order; tuples member by member, including both `Action` operands |
| R010 | Binding declarations; non-expression `Base` objects | Do not visit |
| R011 | Domains, maps, spaces, indices and labels | May inspect, but do not traverse or modify stored expressions |
| R012 | Expressions exposed by a rule | Visit as separate roots; includes newly built map parameterisations |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R013 | Prepared geometry | Geometry is *prepared* when it satisfies core expressions, actions expanded, derivatives expanded and mesh integrals expanded. Reference recipes require the supplier to prepare geometry before exposing it. Constructors do not require this; another strategy may process unprepared geometry with suitable rules. |

Normal-form checks visit these same nodes. Stored dependencies are checked
separately when needed; labels need no lowering. Substitution follows more
dependencies, as described below.

### Substitution and binding

Substitution must respect local names. For example, replacing a free `i` must
not change an `i` bound by an index sum.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R014 | Value operands and labels | Value operands are `AbstractExpression` construction arguments, including those in tuples. Terminals have none. Labels identify objects and follow [Derived indices](expressions.md#derived-indices). Only these nodes bind names. |

| Rule | Binding node | Value operands | Bound names | Bound in |
| ---- | ------------ | -------------- | ----------- | -------- |
| R015 | [`ComponentTensor(e, indices)`, `IndexSum(e, indices)`](expressions.md#componenttensor-and-indexsum) | `e` | The indices of `indices` | `e` and `indices` |
| R016 | [`Form(e, arguments)`](forms.md#form) | `e` | The functions of `arguments` | `e` and `arguments` |
| R017 | [`DirectionalDerivative(f, x, a, h)`](geometry.md#directionalderivative) | `f`, `a`, `h` | $x$ | $f$ and $x$ |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R018 | Substitution traversal | Substitution follows construction dependencies within these scopes and label rules, visiting tuples member by member. Each pass states which operands it transforms. |
| R019 | Fresh bound names | Making a binding node *fresh* for expressions $s_1, \dots, s_n$ means renaming each bound name that occurs free in any $s_k$. Rename it only within its binding scope. This prevents an inserted free name from becoming bound accidentally. |

Substitution replaces $r_1, \dots, r_n$ by $s_1, \dots, s_n$
simultaneously. Within a binding scope, follow these steps.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R020 | Within a binding scope, step 1 | Drop each pair whose key $r_k$ has a bound name free in it. |
| R021 | Within a binding scope, step 2 | Make the binding node fresh for the remaining replacement values $s_k$. |
| R022 | Within a binding scope, step 3 | Substitute the remaining pairs into the operands in that scope. |
| R023 | Outside the binding scope | Operands outside the binding scope use the original pairs. |

To make names fresh for $s_1, \ldots, s_n$, let `body` be the tuple of operands
in the binder's scope. Use these names.

| Rule | Bound name | Fresh name |
| ---- | ---------- | ---------- |
| R024 | Index `i` | `Index(("uflx.rename", i, body, (s_1, ..., s_n)))` |
| R025 | Function `v` in $V$ | `Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))` |

The [derived labels](expressions.md#derived-indices) keep these names distinct
and absent from `body` and the inserted expressions.

`replace(e, mapping)` substitutes each key of the dictionary `mapping` with
its value, simultaneously.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R026 | Matching | Match whole original nodes before their dependencies, respecting binding scopes and label rules. The outermost eligible match wins, even when keys overlap. Do not search inserted values or rebuilt nodes. |
| R027 | Replacement order | Order replacement pairs by the structural key order of [`Base`](rules.md#the-base-object-base), preserving that order through scope filtering and in fresh-name labels. Dictionary insertion order has no effect. |

```python
replace(a + b, {a + b: c, a: d})  # c
replace(a + b, {a: d, d + b: c})  # d + b; no match after rebuilding
replace(a, {a: b, b: c})  # b; inserted values are not searched
```

In these examples, `j` is a fresh index.

| Substitution | Result |
| ------------ | ------ |
| `replace(index_sum(a[i] * c, (i,)), {c: b[i]})` | `index_sum(a[j] * b[i], (j,))`; `i` remains free |
| `replace(form(a, (v,)), {v: w})` | The original form; `v` is bound |
| `replace(DirectionalDerivative(f, u, u, h), {u: u0})` | `DirectionalDerivative(f, u, u0, h)`; the point is outside the binding scope |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R028 | Dependencies and compositions | Substitution follows all construction dependencies, including domains, maps, spaces and derived labels. It commutes with a fixed composition when the static-type choice is unchanged. For example, replacing `A` by `B` in `transpose(A)` gives `transpose(B)`, with the derived indices of `B`. |

### Grammars

A grammar is a set of allowed expressions. A contract check tests membership
of that set (see [New questions](extensions.md#new-questions)).

A grammar pattern inspects a node and objects within a bounded number of links,
without running a pass. It may classify a domain or map without adding that
dependency to the traversal. For a grammar $N$ and patterns $P$, $N \setminus P$
contains the expressions of $N$ with no visited node matching $P$. Combining
contracts means taking the intersection of their grammars.

A recipe accepts expressions in its input grammar and, under its stated
conditions, produces expressions in its output grammar. A rule guarantees only
its own replacement and what it preserves. Extensions must define grammars
and rules for their nodes before those nodes can be lowered to core expressions.

The form layer contains every `Action` and every node with nonempty argument
spaces. The integrand layer contains all other nodes. A contract may restrict
these separately. A `Function` can satisfy core normal form even if its stored
parameterisation contains a reducible gradient. That parameterisation is checked
separately when exposed.

### Representation contracts

| Contract | Required representation |
| -------- | ----------------------- |
| Core expressions | Structurally valid concrete core nodes and downstream terminals |
| Actions expanded | Every `Action` applies a `Function`, downstream terminal, `Gradient`, `PullBack` or another retained `Action` |
| Derivatives expanded | No `DirectionalDerivative` |
| Mesh integrals expanded | No `Integral` over a `Mesh` |
| Gradients expanded | Every remaining `Gradient`, `Divergence` and `Curl` is allowed by [The core normal form](#the-core-normal-form) |
| Pull-backs expanded | Every remaining `PullBack` is allowed by [The core normal form](#the-core-normal-form) |
| Reference integration | Derivatives expanded and mesh integrals expanded |
| Core normal form | The allowed patterns in [The core normal form](#the-core-normal-form) |

A check reports which nodes fail a contract and why, without changing the
expression. Preservation requirements cover domains and operators. Later steps,
including quadrature and tabulation, have their own requirements.

## Actions and derivatives

### Commuting rules

Differentiating a sum gives a sum of derivatives. The table lists this and
other commuting rules. They apply when all other operands are independent of
the operation's target. Specific rules take precedence. Rename bound names
where needed to avoid capturing names in inserted
operands (see [Substitution and binding](#substitution-and-binding)).

| Rule | Node class | Operands | Operations that commute |
| ---- | ---------- | -------- | ----------------------- |
| R029 | `Sum`, `Tensor` | Every operand | `Action`, `DirectionalDerivative`, spatial partial derivative |
| R030 | `Indexed`, `ComponentTensor`, `IndexSum`, `Form`, `DirectionalDerivative` | The expression | `Action`, `DirectionalDerivative`, spatial partial derivative |
| R031 | `Division` | The numerator | `Action`, `DirectionalDerivative`, spatial partial derivative |
| R032 | `Gradient`, `Divergence`, `Curl`, `PushForward`, `PullBack`, `Integral` | The expression | `Action`, `DirectionalDerivative` |
| R033 | `Conjugate`, `RealPart`, `ImaginaryPart` | The expression | `DirectionalDerivative`, spatial partial derivative |

Directional differentiation uses a real parameter and holds domains fixed.
`Action` has its own [conjugation rule](#applying-actions); `RealPart` and
`ImaginaryPart` have no argument spaces. Actions on gradients and pull-backs
stay outside them, as specified in their rule tables.

### Applying actions

An action supplies a value for a form's argument. This pass moves actions
into operands and substitutes values into explicit forms.

| `apply_actions(f)` |  |
| ------------------ | - |
| Evaluates          | `Action` |
| Input contract     | Core expressions |
| Output contract    | Actions expanded |
| Strategy | Inner actions first; first matching row; revisit created actions, including those exposed by substitution |

Downstream terminals may declare argument spaces, for example for assembled
vectors or matrices; their libraries evaluate the retained actions.

| Rule | Form $f$ | Rewrite of `Action(f, p, g)` |
| ---- | -------- | ---------------------------- |
| R034 | `Form(e, arguments)`, $p$ at least the arity $n$ of $e$ | `Form(e, others)` with $g$ substituted for the argument at $p - n$, or $e$ with $g$ substituted if no entries remain |
| R035 | `Conjugate(e)` | `conjugate(action(e, p, conjugate(g)))`, with conjugation member by member for a `tuple` $g$ |
| R036 | `Product(a, b)`, with $n$ the arity of $a$ | `Product(Action(a, p, g), b)` for $p < n$, otherwise `Product(a, Action(b, p - n, g))` |
| R037 | `Gradient`, `PullBack` | Remains; the operation transforms the form's output |
| R038 | A commuting node class (see [Commuting rules](#commuting-rules)) | The node with the action moved into the listed operands |
| R039 | `Zero` | Zero with the type of the action |
| R040 | `Function` or a downstream terminal | Remains, for the supplying library |
| R041 | `Action` | Remains, as a further application of a remaining action |
| R042 | Any other class of the core | Cannot occur, because the argument spaces are `()` |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R043 | Binding during action substitution | Substitution preserves all enclosing and inserted binders. In particular, commuting through `DirectionalDerivative` makes its variable fresh for `g`; occurrences in `g` are not differentiated. |

### Derivatives of the node classes

Directional and spatial derivatives use this table after their specific rules.
Here $\delta b$ is the derivative of $b$. Each pass defines independence from
its differentiation target. Extensions can add rules here.

| Rule | Expression | Derivative |
| ---- | ---------- | ---------- |
| R044 | `Product(b, c)` | $\delta b \, c + b \, \delta c$ |
| R045 | `Division(b, c)` | $\delta b / c - (b \, \delta c) / (c c)$ |
| R046 | `Power(b, c)`, $c$ constant | $c \, b^{c - 1} \, \delta b$ |
| R047 | `Power(b, c)` | $b^c (\log(b) \, \delta c + c \, \delta b / b)$ |
| R048 | `Absolute(b)` | $\operatorname{Re}(\bar b \, \delta b) / \lvert b \rvert$, with `real_part`, so the derivative is real |
| R049 | `SquareRoot(b)` | $\delta b / (2 \sqrt b)$ |
| R050 | `Logarithm(b)` | $\delta b / b$ |
| R051 | `Determinant(A)`, $A$ invertible | $\det(A)\sum_{i,j}(A^{-1})_{ij}\,\delta A_{ji}$, with `inverse`, `index_sum` and derived indices $i,j$ |
| R052 | `Inverse(A)`, equal axis extents, including symbolic extents | `-inverse(A) @ delta_A @ inverse(A)`, with `delta_A` = $\delta A$ |
| R053 | `Inverse(A)`, rectangular shape | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |
| R054 | `Conjugate(b)`, `RealPart(b)`, `ImaginaryPart(b)` | The same node applied to $\delta b$ |
| R055 | A commuting node class (see [Commuting rules](#commuting-rules)) | The node with the derivative moved into the listed operands |

The constant-exponent power rule avoids a logarithm, permitting negative real
bases where differentiable. Results are core expressions.

### Applying derivatives

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` and actions created by its evaluation |
| Input contract         | Core expressions; actions expanded |
| Output contract        | Actions and derivatives expanded |
| Strategy | Inner derivatives first; first matching row; alternate with `apply_actions` until both output contracts hold |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R056 | Derivative substitution | Write $\mathrm{D}b$ for `DirectionalDerivative(b, x, a, h)`. Match original operands before substituting `[x := a]` into undifferentiated operands; insert `h` unchanged. |

| Rule | Expression $f$ | Rewrite of the derivative at $a$ |
| ---- | -------------- | -------------------------------- |
| R057 | No function of $x$ free in $f$ | Zero with the type of the derivative |
| R058 | $f$ equal to the member $x_k$ of $x$ | $h_k$, or $h$ for a function $x$ |
| R059 | `Action(b, p, g)` | `Action(Db, p, g) + Action(b, p, Dg)`, with undifferentiated operands evaluated at $a$ |
| R060 | `DirectionalDerivative` | Defer until the inner derivative is expanded |
| R061 | Any other expression | The row of [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \mathrm{D} b$ |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| R062 | Dependence on the derivative variable | To test whether an expression depends on `x`, follow value operands and binding scopes, holding domains and maps fixed. A `Function` contributes itself, not names in its label or space. A form binding `x` has no free `x`. Tuple arguments vary member by member, including conjugated entries. Dual directions and substitution points may be symbolic forms whose actions the recipe evaluates. |

For real scalar functions `u` and `a` in a space of constants,
`directional_derivative(u * u, u, a, u)` rewrites to `u * a + a * u`.

For a real scalar function `u` in a space of constants and a constant array
`a` of shape `(m,)`, `directional_derivative(u, u, a[i], Integer(1))` lowers
to `Integer(1) + Zero((), ((i, m),), (), real)`. Its value is one for every `i`.
Type completion retains indices even from unused members of tuple points
and directions.

## Spatial lowering

### Expanding integrals

Mesh integration becomes a sum over cells, with each cell integral expressed
on its reference cell. The density accounts for the change of coordinates.

| `expand_integrals(f)` |  |
| --------------------- | - |
| Evaluates             | `Integral` over a mesh |
| Input contract        | Core expressions |
| Output contract       | Mesh integrals expanded |
| Strategy | Inner mesh integrals first |

For group $j$ of `D.topology`, let $K_j$ be the reference cell and $c_j$ =
`cell_index(D.topology, j)`.

To choose the internal index without capturing $f$, form the template
`index_sum(pull_back(Integer(1), cell_map(D, j)), (c_j,))`, make its binder
fresh for $f$ using [Substitution and binding](#substitution-and-binding), and
let $e_j$ be its resulting bound index. This retains $c_j$ when no renaming is
needed, preserving shared geometry for ordinary integrands. Set
`D_j = cell_map(D, j, e_j)` and
`piece_j = domain_map(K_j, D_j.parameterisation, D.ambient_domain)`.

Let `K` contain one derived `Index` for each value axis of $f$, with `()`
for a scalar, and set `f_K = indexed(f, K)`.

| Rule | Integral | Rewrite |
| ---- | -------- | ------- |
| R063 | `Integral(f, D)` | `component_tensor(sum_j, K)`, where `sum_j` is the sum over groups $j$ of `index_sum(integral(g_j * density(D_j), K_j), (e_j,))` |

$g_j$ is `pull_back(f_K, piece_j)` when $f$ has domain `D.ambient_domain`,
and `pull_back(f_K, D_j)` otherwise, including domain-free $f$.
The pull-back introduces the internal cell index; `IndexSum` sums over it and
retains every free index of $f$ and the component indices `K`. The outer
`component_tensor` restores the value shape, and is the identity for scalar
$f$. The remaining integrals have coordinate domains.

Expanding the same input again gives the same graph. Changing the order of
substitution and expansion may change bound names when renaming is needed.
The value and free indices must stay the same.

### Applying gradients

This pass applies spatial calculus rules while keeping derivatives that need
the supplying library or a later change of coordinates.

| `apply_gradients(f)` |  |
| -------------------- | - |
| Evaluates            | `Gradient`, `Divergence` and `Curl` |
| Input contract       | Core expressions; actions and derivatives expanded |
| Output contract      | Actions, derivatives and gradients expanded |
| Strategy | Inner expressions first; first matching row; evaluate exposed actions with `apply_actions` and revisit affected expressions |

A cell is a `CellMap` of a topology or mesh. Rows marked “Boundary” keep an
expression for later; the divergence and curl rows preserve the Piola identities.
The family column groups rules for selection by a strategy. Families are not
separate public passes.

$\partial_j b$ denotes component $j$ of the gradient of $b$ on $P$, with
$b$'s value shape and the additional free index $j$ of ambient extent $n$.
Gradient reconstruction uses derived component indices $A$ and $j$, one
index in $A$ per value axis of $f$. Empty lists use the factory identities.

| Rule | Family | Expression or operation | Rewrite |
| ---- | ------ | ----------------------- | ------- |
| R064 | Action transport | `Gradient(Action(B, p, g), P)` | `Action(Gradient(B, P), p, g)` |
| R065 | Boundary | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)` | Remains for the Piola rule of a cell pull-back |
| R066 | Boundary | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$ | Remains for the Piola rule of a cell pull-back |
| R067 | Differential expansion | Other `Divergence` or `Curl` | Its composition of gradients from [Divergence and Curl](geometry.md#divergence-and-curl) |
| R068 | Gradient reconstruction | Other `Gradient(f, P)` | `component_tensor(∂_j f[A], A + (j,))`, using the component rules below |
| R069 | Gradient rules | $\partial_j b$, $b$ with no domain | Zero with $b$'s shape and argument spaces, and its free indices followed by $j$ |
| R070 | Gradient rules | $\partial_j$ `Coordinate(P)`, $P$ of its ambient dimension | `component_tensor(I[i, j], (i,))`, with $I$ the identity of shape `(n, n)` and a derived $i$ |
| R071 | Boundary | $\partial_j b$, $b$ a terminal, `Gradient`, `PushForward`, Piola-matched node, or `PullBack` by a cell or `SideMap` | `component_tensor(gradient(b, P)[B + (j,)], B)`, with one derived index in `B` per value axis of $b$; the gradient remains |
| R072 | Gradient rules | $\partial_j$ `PullBack(g, φ)`, φ a `DomainMap` | `component_tensor(index_sum(H[A + (l,)] * J[l, j], (l,)), A)`, with `H = pull_back(gradient(g, Q), φ)`, `J = jacobian(φ)`, and $Q$ the ambient domain of φ |
| R073 | Action transport | $\partial_j$ `Action(B, p, g)` | `Action(∂_j B, p, g)`, with $g$ unchanged |
| R074 | Gradient rules | $\partial_j b$, any other $b$ | [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \partial_j b$ |

In the chain rule, `A` indexes the value axes of `g` and `l` is a distinct
contraction index; scalar `g` uses `A == ()`. On embedded meshes the same rules
apply to the tangential gradient. The coordinate gradient remains a projection
until the pull-back to a cell. The supplying library evaluates remaining terminal
gradients. The next call to `apply_gradients` handles gradients exposed by pull-backs.

### Expanding pull-backs

Pull-backs express values in a map's parameter coordinates. This pass moves
them through arithmetic and applies the geometry and mapping rules below.
Boundary rows keep reads that the supplying library must evaluate.

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

| Rule | Family | `PullBack(f, ψ)` | Rewrite |
| ---- | ------ | ---------------- | ------- |
| R075 | Action transport | `Action(A, p, g)` | `Action(PullBack(A, ψ), p, g)`, with $g$ unchanged |
| R076 | Binding | `Form(e, arguments)` | Make the form fresh for $\psi$, obtaining `Form(e_new, arguments_new)`; return `Form(PullBack(e_new, ψ), arguments_new)`, preserving argument spaces |
| R077 | Boundary | $\psi$ a `SideMap`, $f$ with a domain | Remains until a pull-back by a cell of $P$ |
| R078 | Pointwise | $f$ with no domain | $f$, with [type completion](#type-compatibility-of-rewrites) to retain the free indices of $\psi$ |
| R079 | Geometry | `Coordinate(Q)` | The parameterisation of $\psi$ |
| R080 | Restriction | `PullBack(g, SideMap(X, Y, s))`, $\psi$ = `cell_map(Y, i, q)` | `PullBack(PullBack(g, B_selected), A_selected)` |
| R081 | Boundary | `PullBack(g, φ)`, φ a `SideMap` or a cell | Remains, so a pull-back nests only inside a read on a cell |
| R082 | Map composition | `PullBack(g, φ)`, φ a `DomainMap` | `PullBack(g, domain_map(P, pull_back(φ.parameterisation, ψ), φ.ambient_domain))`, the pull-back by the image of $\psi$ under φ (see [`PullBack`](geometry.md#pullback)) |
| R083 | Mapped values | `PushForward(f, D, m)`, $\psi$ a cell of $D$ | $m_\psi(\hat f)$ |
| R084 | Piola | `Divergence(PushForward(f, D, m), D)` of the identity of Piola, $\psi$ a cell of $D$ | $m'_\psi$ of `divergence(f̂, K_j)`, with $m'$ = `transformation_law((), 1)` |
| R085 | Piola | `Curl(PushForward(f, D, m), D)` of the identity of Piola, $\psi$ a cell of $D$ | $m'_\psi$ of `curl(f̂, K_j)`, with $m'$ = `transformation_law(("contravariant",), 1)` for $n = 3$ and `transformation_law((), 1)` for $n = 2$ |
| R086 | Boundary | `PushForward`, `Divergence` or `Curl` | Remains, read at the points of $P$ |
| R087 | Gradient transport | `Gradient(g, Q)`, $\psi$ a cell of $Q$ | $m_\psi(\nabla_P(\psi^* g))$, with $m$ = `transformation_law(("covariant",), 0)` |
| R088 | Boundary | `Gradient(g, Q)` | Remains, the gradient on $Q$ read at the points of $P$ |
| R089 | Geometry | `Normal(SideMap(D, E, s))`, $\psi$ = `cell_map(E, i, q)` | $w / \lvert w \rvert$ of [`Normal`](geometry.md#normal) |
| R090 | Boundary | A terminal | Remains, for the library of the terminal |
| R091 | Pointwise | Any other node class | Make its binders fresh for $\psi$; rebuild from the pull-backs of its expression operands |

Map composition processes its new parameterisation as a separate root under
the [geometry preparation contract](#evaluation-traversal). Cell maps retain
explicit selection; unselected side restrictions remain. Piola rows read no
derivative of $J$. Only cell transport inverts a Jacobian; other maps use the
[gradient chain rule](#applying-gradients). Libraries evaluate remaining cell
reads under [Functions in a space](forms.md#functions-in-a-space).

## The core normal form

The core normal form combines
the core expressions, actions expanded, derivatives expanded, mesh integrals
expanded, gradients expanded and pull-backs expanded contracts. The table lists
the nodes allowed after `lower`.

| Chapter | Node class | Allowed in the core normal form |
| ------- | ---------- | ------------------------------- |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Integer`, `Float`, `Complex`, `Zero` | Always |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Sum`, `Product`, `Division`, `Power`, `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm` | Always |
| [Tensors and indexing](expressions.md#tensors-and-indexing) | `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | Always; index simplification is left to the compiler |
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

The normal-form check follows the [evaluation traversal](#evaluation-traversal);
it does not check stored geometry or labels.

Assembly code receives the form layer (see [Grammars](#grammars)), including
sums and products of forms, real and complex zeros, nested forms, and actions
on supplied or spatially transformed forms. The static type gives the argument
order (see [`Form`](forms.md#form)).

A form compiler receives the integrand layer, consisting of integrals over reference cells
and points, with arithmetic and index notation over constants, reference
coordinates and the library reads of [Functions in a space](forms.md#functions-in-a-space).

```{dropdown} Reusing geometry calculations
Jacobians, inverse Jacobians, densities and side normals built by these passes
depend only on domain maps and remain subexpressions of the result. A compiler
can compute them once for given map data, such as mesh coordinates, and reuse
them in every integral that needs them.
```

## Standard lowering recipe

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `lower(f)` | | `f: AbstractExpression` | An equivalent expression in core normal form |

| `lower` condition | Requirement or result |
| ----------------- | --------------------- |
| Input | Core expressions and [prepared geometry](#evaluation-traversal); extension nodes need an extension recipe first |
| Unsupported rewrite | `TypeError` |
| Output violates core normal form | `ValueError`, identifying offending nodes |

It runs these passes in order.

| Recipe step | Required result |
| ----------- | --------------- |
| [`apply_actions`](#applying-actions) | Actions expanded |
| [`apply_derivatives`](#applying-derivatives), including action cleanup | Actions and derivatives expanded |
| [`expand_integrals`](#expanding-integrals) | Also mesh integrals expanded |
| Repeat [`apply_gradients`](#applying-gradients), then [`expand_pull_backs`](#expanding-pull-backs), with action cleanup, until neither changes the expression | Core normal form |

This order evaluates gradients on their stated domains before changing coordinates.
Full-dimensional `gradient(coordinate(D), D)` becomes the identity before
introducing Jacobians. Expressions matching Piola identities remain until a
matching cell pull-back.

Other recipes may choose individual rules or rule families, provided their
conditions hold and geometry is prepared as required. For example,
[shape differentiation](extensions.md#shape-derivatives) lowers spatially before
and after differentiation. Index-notation simplification, such as reducing an
`Indexed` of a `ComponentTensor`, remains a compiler operation.
