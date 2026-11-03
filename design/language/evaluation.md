# Evaluation and lowering

## Evaluation

A lazy node class states its value in its section. A pass may rewrite it as
an expression of other node classes with that value. The lazy classes are
`Action`, `DirectionalDerivative`, `Integral` over a mesh, `PullBack`,
`Gradient`, `Divergence`, `Curl` and `PushForward`.

A rewrite preserves value under the [mathematical admissibility](rules.md#mathematical-admissibility)
convention and preserves [type compatibility](#type-compatibility-of-rewrites),
using the renaming of [Substitution and binding](#substitution-and-binding).
Value preservation is mathematical; it does not guarantee identical
floating-point results after rewriting.

The evaluation rules do not require algebraic simplification; zero terms may
remain in the graph. Implementations may apply value- and type-preserving
optimizations under the mathematical admissibility convention. No rule makes
zero times an undefined value defined.

`Action` and `DirectionalDerivative` commute with the following node classes
in the listed operands, when the other operands do not depend on what the
evaluated node changes:

| Node class | Operands |
| ---------- | -------- |
| `Sum`, `Tensor` | Every operand |
| `Indexed`, `ComponentTensor`, `IndexSum`, `Form`, `DirectionalDerivative` | The expression |
| `Division` | The numerator |
| `Gradient`, `Divergence`, `Curl`, `PushForward`, `PullBack`, `Integral` | The expression |

`DirectionalDerivative` also commutes with `Conjugate`, `RealPart` and
`ImaginaryPart`, because its parameter is real. `Action` has a separate rule
for `Conjugate` (see [Applying actions](#applying-actions)); `RealPart` and
`ImaginaryPart` have no argument spaces.

The value of a commuting node is linear in the listed operands, and the
evaluated node acts on the value of an operand alone: an action substitutes into
a bound field or functional, and a directional derivative
differentiates in a real parameter with the domains fixed. Moving an evaluated
node into a binding node first makes the binding node fresh for the other
operands of the evaluated node (see [Substitution and
binding](#substitution-and-binding)). The partial derivative of
[`apply_gradients`](#applying-gradients) commutes with the same node classes
except `Gradient`, `Divergence`, `Curl`, `PushForward`, `PullBack` and
`Integral`, which the rows of that pass match first; it also commutes with
`Conjugate`, `RealPart` and `ImaginaryPart`.

For `Action`, the `Gradient` and `PullBack` identities are evaluated toward
actions on transformed forms, as in [`Action`](forms.md#action).

### Transformations, recipes and contracts

A transformation changes an expression using the rules below. A recipe
coordinates transformations. A contract describes an accepted representation
and can be checked without transforming it. These are separate responsibilities;
a contract does not choose a recipe.

Consumers may compose transformations and stop at a suitable boundary:

| Contract | Required representation |
| -------- | ----------------------- |
| Derivatives expanded | No `DirectionalDerivative` in the expression traversal; actions, mesh integrals and spatial operators may remain |
| Reference integration | Derivatives expanded, with no `Integral` over a `Mesh`; other coordinate-domain integrals and spatial operators may remain |
| Core normal form | The node classes and patterns of [The core normal form](#the-core-normal-form) |

A check reports offending nodes and reasons, or no violations. It assumes valid
nodes and checks only the [evaluation traversal](#evaluation-traversal), not
mathematical declarations or later compiler capabilities. Quadrature, tabulation
and geometry processing impose their own requirements at their own boundaries.
No target class or checking API is prescribed.

An output contract does not establish that a transformation preserved particular
input structure. A symbolic consumer that needs the original domains and spatial
operators retained must also require that preservation from its chosen recipe.

### Standard lowering recipe

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `lower(f)` | | `f: AbstractExpression` | An equivalent expression in core normal form |

`lower(f)` runs the standard recipe and checks its output. It raises `TypeError`
for a non-expression input or an unsupported rewrite, and `ValueError` if the
result violates core normal form, identifying the offending nodes. It assumes
the geometry precondition of [Evaluation traversal](#evaluation-traversal).
Extension nodes require an extension recipe before entering this core operation.

The standard recipe takes a core expression to core normal form. The table gives
one schedule and names its intermediate grammars for the rule tables below.
Consumers need not follow this schedule when targeting another contract, but must
satisfy the preconditions of the transformations they use.

| Pass | Evaluates | Output grammar |
| ---- | --------- | -------------- |
| [`apply_actions`](#applying-actions) | `Action`, except actions on a `Function`, downstream terminal, `Gradient`, `PullBack` or another remaining `Action` | $N_1$ |
| [`apply_derivatives`](#applying-derivatives) | `DirectionalDerivative`, then any actions it creates | $N_2$ |
| [`expand_integrals`](#expanding-integrals) | `Integral` over a mesh | $N_3$ |
| The lowering stage | `Gradient`, `Divergence`, `Curl`, `PushForward` and `PullBack` | $N_4$, [the core normal form](#the-core-normal-form) |

In this recipe, the first three passes keep this order. `apply_actions` evaluates actions
created by substitution, and `apply_derivatives` uses it to evaluate actions
created by differentiation or substitution at the derivative's point.
The lowering stage repeats [`apply_gradients`](#applying-gradients)
and then [`expand_pull_backs`](#expanding-pull-backs) until neither rewrites.
Both passes evaluate actions they create or expose with `apply_actions` and
revisit affected expressions with their own rules.
`apply_gradients` runs first, so a gradient is evaluated on the domain on
which the form states the gradient, such as `gradient(coordinate(D), D)`, the
identity matrix for $D$ of the dimension of its ambient coordinate domain, and
not as the product of $J$ and $K$ on the reference cell.

The recipe checks its output against core normal form. Constructs allowed there
remain for later processing; a construct outside that contract is reported,
not silently accepted as a completed lowering. Other recipes may have weaker
output contracts. A transformation's own postcondition is the one in its rule
table, not necessarily the contract of the complete recipe.

The passes use the [evaluation traversal](#evaluation-traversal). A substitution
of indices also enters domains, domain maps and spaces. Extensions coordinate
their transformations with core transformations through their own recipes. For
example, [shape differentiation](extensions.md#shape-derivatives) needs geometry exposed by
spatial lowering before differentiating, and spatial lowering again afterwards.

A rewrite of index notation that keeps the value, such as the reduction of an
`Indexed` of a `ComponentTensor`, is a pass of the compiler, not of the core.

## Rewriting and substitution

### Type compatibility of rewrites

A replacement $r$ is type-compatible with the expression $e$ it replaces if:

- Their value shapes, scalar types and ordered argument spaces are equal.
- Their free indices have the same names and extents, up to order and the
  consistent renaming of bound names.
- Their domains are equal, or $r$ has no domain. In the latter case,
  $r$ represents a constant field on the domain of $e$, for each assignment
  of the free indices and form arguments.

The replacement must preserve the value on the original domain. A pull-back
of `Integer(2)` may become `Integer(2)`, and the gradient of coordinates on a
full-dimensional domain may become a literal identity matrix. These results
lose domain; a consumer needing the original domain must retain it.

Before returning a table's result `r` or rebuilding its parent, complete its
type against the replaced expression `e`. Let

```python
z = Zero(e.value_shape, e.free_indices, e.argument_spaces, e.scalar_type)
```

Apply these steps in order, using tensor arithmetic:

- If `e` is complex and `r` is real, promote `r` to
  `Complex(Integer(1), Integer(0)) * r`.
- If `r` lacks any free index of `e`, replace `r` by `r + z`.

The raw result must already have the required shape and argument spaces,
matching extents for retained indices and no extra free indices. A real `e`
requires a real result. A row giving zero uses `z`; this is what the tables
mean by zero with the type of an expression.

Evaluation of `Action` and `DirectionalDerivative` also promotes real values
before substituting them for complex functions, member by member for tuples.
General structural substitution with `replace` does not insert this promotion.
For example, substituting a real negative point into `square_root(u)` for
complex `u` must preserve complex arithmetic. Promotion after rebuilding the
square root would be too late. Likewise, missing family indices must be
restored before rebuilding an enclosing `ComponentTensor` or `IndexSum`.

### Evaluation traversal

Evaluation visits the root expression and recursively visits its value
operands in construction order, as defined in [Substitution and
binding](#substitution-and-binding). Traversal stops at non-expression `Base`
objects and does not visit binding declarations.

Thus evaluation visits both expression arguments of an `Action`, including
a tuple argument, but does not enter domains, maps, spaces, indices or their
labels. Rules may inspect these dependencies and expose expressions from
them: a pull-back of coordinates returns a map's parameterisation, for
example. Exposed expressions are visited; stored geometry is unchanged.

Each expression exposed from geometry must belong to $N_3$ as a separate
root. The supplier runs the first three passes before building the geometry
if necessary. This is a pass precondition, not a constructor restriction;
it also applies to parameterisations of maps that rewrites build. The lowering
stage evaluates exposed gradients and pull-backs.

The enclosing expression's normal form makes no claim about expressions
retained in geometry, spaces or metadata. Consumers evaluate these separately
when needed. Labels have no lowering-grammar requirement. Substitution still
follows its dependency and binding rules, including rebuilding derived labels
without evaluating their contents.

### Substitution and binding

The following rules define value operands, binding declarations and scopes
for evaluation, substitution and differentiation:

- Terminals have no value operands.
- Except for the nodes below, value operands are the `AbstractExpression`
  construction arguments, including those nested in `tuple`s; no names are bound.
- Labels are identity metadata, not value operands. Substitution follows the
  reconstruction rules of [Derived indices](expressions.md#derived-indices).

| Binding node | Value operands | Bound names | Bound in |
| ------------ | -------------- | ----------- | -------- |
| [`ComponentTensor(e, indices)`, `IndexSum(e, indices)`](expressions.md#componenttensor-and-indexsum) | `e` | The indices of `indices` | `e` and `indices` |
| [`Form(e, arguments)`](forms.md#form) | `e` | The functions of `arguments` | `e` and `arguments` |
| [`DirectionalDerivative(f, x, a, h)`](expressions.md#directionalderivative) | `f`, `a`, `h` | $x$ | $f$ and $x$ |

Tuple-valued operands are visited member by member. Substitution follows
construction dependencies subject to these scopes and the label rules;
individual passes specify which value operands they transform.

A binding node is made fresh for expressions $s_1, \dots, s_n$ by renaming,
in the operands in which the names are bound, each bound name that is free in
some $s_k$, so that the binding node does not capture the name.

A substitution replaces expressions $r_1, \dots, r_n$ by $s_1, \dots, s_n$
simultaneously, so a substituted $s_k$ is not searched for the other $r_l$.
Inside a binding node, the substitution goes into the other operands as
anywhere else. For the operands in which the names are bound, the substitution
drops each pair whose $r_k$ has a bound name free in it, because every
occurrence of that $r_k$ there refers to the bound name, makes the binding node
fresh for the remaining $s_k$, and substitutes the remaining pairs.

When a binding node is made fresh for $s_1, \dots, s_n$, the new name of a
bound index `i` is the derived index
`Index(("uflx.rename", i, body, (s_1, ..., s_n)))`, and the new
name of a bound function `v` in a space $V$ is
`Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))`, with `body` the
`tuple` of the operands in which the names are bound. The label refers to
`body` and to each $s_k$, so the new name occurs in neither (see [Derived
indices](expressions.md#derived-indices)). The bound names of a binding node are distinct,
so the new names are distinct.

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

For example, substituting an expression with a free index requires renaming:

```python
e = index_sum(a[i] * c, (i,))
result = replace(e, {c: b[i]})
```

With `j` the renamed index, the result is:

```python
index_sum(a[j] * b[i], (j,))
```

The index `i` remains free. Without renaming, it would be captured by the sum.

A bound function is not replaced:

```python
F = form(a, (v,))
replace(F, {v: w})  # Leaves the bound v unchanged.
```

The point of a derivative is not bound, even when it is the same function as
the variable:

```python
d = DirectionalDerivative(f, u, u, h)
result = replace(d, {u: u0})
```

The result is the derivative at `u0`:

```python
DirectionalDerivative(f, u, u0, h)
```

### Grammars

A grammar is a set of expressions and specifies a structural contract. A check
asks whether an expression belongs to that set (see [New questions](extensions.md#new-questions)).
The numbered grammars below describe the standard recipe; they are not a sequence
that every consumer must traverse.

- For the lowering grammars, an expression node of an expression $e$ is a
  node visited by the [evaluation traversal](#evaluation-traversal), including
  $e$. The grammar $N_0$ is the set of the
  expressions $e$ in which every expression node is an instance of a concrete
  class of the core or a terminal of a downstream library (see [Leaves and
  subclassing](rules.md#leaves-and-subclassing)). Every object is structurally valid,
  so $N_0$ needs no other structural condition. A node class of a language
  extension is not in $N_0$, and an extension that adds node classes states
  the grammars that contain them.
- A pass maps a source grammar into a target grammar if, for every
  expression of the source grammar satisfying the pass's preconditions, the
  result belongs to the target grammar. The core pass contracts also require
  the geometry precondition of [Evaluation traversal](#evaluation-traversal).
- A grammar pattern is a condition on one expression node that reads only
  the node and the objects within a bounded number of edges. It never reads
  the result of a pass. It may inspect a domain or map to classify that node;
  this does not add the expressions in that dependency to the traversal.
- For a grammar $N$ and a set of grammar patterns $P$, the delta $N \setminus P$ is
  the set of the expressions of $N$ with no expression node that matches a
  pattern of $P$. An intermediate grammar is written as a delta of an earlier
  grammar.
- The core normal form is the output contract of the standard recipe, written
  as a list of the node classes and the grammar patterns that it allows. A
  form compiler receives the integrand layer of an expression in the core
  normal form and continues with passes of its own, and the code that
  assembles receives the form layer.

A grammar may restrict two layers of an expression separately. The form layer
consists of every `Action` and every node with argument spaces other than
`()`. The integrand layer consists of the other nodes.

Each pass states its input and output grammars. These are conditional guarantees:
for example, `expand_integrals` removes mesh integrals, but reaches the reference
integration contract only when no directional derivatives remain. A consumer
checks the contract it needs rather than inferring success from a pass's name.

For example, a `Function` is in $N_4$ even if its domain's parameterisation
contains a reducible `Gradient`: traversal stops at the function. The
parameterisation is checked separately when exposed or evaluated.

## Actions and derivatives

### Applying actions

| `apply_actions(f)` |  |
| ------------------ | - |
| Evaluates          | `Action` |
| Input grammar      | $N_0$ |
| Output grammar     | $N_1 = N_0 \setminus \{$an `Action` whose form is not a `Function`, downstream terminal, `Gradient`, `PullBack` or `Action`$\}$ |

The pass rewrites every `Action(f, p, g)`, inner actions first, with the first
row that matches the form $f$, and keeps the value of [`Action`](forms.md#action).
It applies the rules again to actions created by a rewrite, including those
created by substituting a symbolic form for a dual function.

A downstream terminal may declare argument spaces, for example to represent
an assembled vector or matrix held by its supplying library. Actions on it
remain for that library to evaluate. The terminal must satisfy the ordinary
argument-space and action contracts.

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

- `DirectionalDerivative` binds its variable, so moving the action into it
  renames the variable if the variable is free in $g$, and the $x$ in $g$ is
  not differentiated.
- A substituted form is a value with its own binders. Substitution preserves
  their scope and the scope of the enclosing form. The pass finishes when
  every remaining action applies a `Function`, a downstream terminal, a
  `Gradient`, a `PullBack` or another remaining action. Actions on `Gradient` and `PullBack` stay
  outside these operations, so spatial lowering does not rewrite in a cycle.

### Derivatives of the node classes

The passes that differentiate share one table, the derivative of each node
class in terms of the derivatives of its operands. In the table, $\delta b$ is
the derivative of the operand $b$ in the pass, and an operand is constant for
the pass if the first row of the pass matches the operand. Each pass matches
its own rows first, for the terminals of the expression, and then the rows below.

| Expression | Derivative |
| ---------- | ---------- |
| `Product(b, c)` | $\delta b \, c + b \, \delta c$ |
| `Division(b, c)` | $\delta b / c - (b \, \delta c) / (c c)$ |
| `Power(b, c)`, $c$ constant | $c \, b^{c - 1} \, \delta b$ |
| `Power(b, c)` | $b^c (\log(b) \, \delta c + c \, \delta b / b)$ |
| `Absolute(b)` | $\operatorname{Re}(\bar b \, \delta b) / \lvert b \rvert$, with `real_part`, so the derivative is real |
| `SquareRoot(b)` | $\delta b / (2 \sqrt b)$ |
| `Logarithm(b)` | $\delta b / b$ |
| `Determinant(A)` | $\det(A)\sum_{i,j}(A^{-1})_{ij}\,\delta A_{ji}$, with `inverse`, `index_sum` and derived indices $i,j$ |
| `Inverse(A)`, square shape | `-inverse(A) @ delta_A @ inverse(A)`, with `delta_A` = $\delta A$ |
| `Inverse(A)`, rectangular shape | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |
| `Conjugate(b)`, `RealPart(b)`, `ImaginaryPart(b)` | The same node applied to $\delta b$ |
| A commuting node class (see [Evaluation](#evaluation)) | The node with the derivative moved into the listed operands |

- The `Determinant` row applies where $A$ is invertible.
- The square `Inverse` row applies when the two value-shape extents are equal,
  including equal symbolic extents; otherwise use the rectangular row.
- The constant-exponent row for `Power` avoids the logarithm of the general
  row, allowing negative real bases where the power is differentiable.
- A commuting row applies only to an operand that no row of the pass
  matches first, so each pass decides through its own rows which node classes
  the derivative moves into, as [`apply_gradients`](#applying-gradients) does
  for `Gradient`, `PushForward` and `PullBack`.
- A language extension that adds a node class adds the derivative of the
  class to this table, and every pass that differentiates then uses the row.
- The derivative of each row is a composition of node classes of the core,
  with actions evaluated as in [Applying actions](#applying-actions).

### Applying derivatives

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` and actions created by its evaluation |
| Input grammar          | $N_1$ |
| Output grammar         | $N_2 = N_1 \setminus \{$a `DirectionalDerivative`$\}$ |

The pass rewrites every `DirectionalDerivative(f, x, a, h)`, inner derivatives
first, with the first row that matches $f$. $\mathrm{D} b$ is short for
`DirectionalDerivative(b, x, a, h)`. The rules give the derivative directly
at $a$: undifferentiated operands in the derivative formulas use the
substitution $[x := a]$ of [`DirectionalDerivative`](expressions.md#directionalderivative),
while $h$ is inserted unchanged. Rules match the original operands, before
substitution. The pass applies `apply_actions` to the result, so the output
remains in $N_1$ and contains no `DirectionalDerivative`.

| Expression $f$ | Rewrite of the derivative at $a$ |
| -------------- | -------------------------------- |
| No function of $x$ free in $f$ | Zero with the type of the derivative |
| $f$ equal to the member $x_k$ of $x$ | $h_k$, or $h$ for a function $x$ |
| `Action(b, p, g)` | `Action(Db, p, g) + Action(b, p, Dg)`, with undifferentiated operands evaluated at $a$ |
| Any other expression | The row of [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \mathrm{D} b$ |
| `DirectionalDerivative` | Cannot occur, because the inner derivative is rewritten first |

- The first row comes first, so no derivative is built for an operand that
  does not depend on $x$, and a terminal other than a function of $x$ needs no
  row. An operand is constant for the pass if no function of $x$ is free in
  the operand.
- This test follows value operands and respects binding scopes. A `Function`
  contributes itself, not occurrences in its label or space. Domains and maps
  are held fixed, as in the definition of `DirectionalDerivative`.
- A `Form` that binds $x$ has no free $x$, so the first row gives zero.
- A dual variable can occur in the functional operand of a pairing. Its
  direction and substitution point may be symbolic forms; their actions
  are evaluated by `apply_actions`.
- The `Action` row varies both expression operands, member by member for a
  `tuple` argument. A conjugated argument varies through `Conjugate`.

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
| Input grammar         | $N_2$ |
| Output grammar        | $N_3 = N_2 \setminus \{$an `Integral` whose domain is a `Mesh`$\}$ |

The pass rewrites mesh integrals, inner integrals first. For group $j$ of
`D.topology`, let $K_j$ be the reference cell and $c_j$ =
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
| Input grammar        | $N_3$, in the lowering stage |
| Output grammar       | $N_3 \setminus G$, with $G$ the grammar patterns of the rows below that rewrite |

A cell is a `CellMap` of a topology or mesh. The identity of Piola
matches `Divergence(PushForward(f, D, m), D)` with $m$ =
`transformation_law(("contravariant",), 1)`, and `Curl(PushForward(f, D, m), D)`
with $m$ = `transformation_law(("covariant",), 0)` and $k = n$ (see
[`PushForward`](geometry.md#pushforward)).

The pass rewrites `Gradient(Action(A, p, g), P)` as
`Action(Gradient(A, P), p, g)`. It rewrites every `Divergence` and `Curl` that
the identity of Piola does not match as its composition of `gradient(u, P)`,
and every other `Gradient(f, P)` that a row below rewrites, inner gradients first, as
`component_tensor(∂_j f[a_1, …, a_r], (a_1, …, a_r, j))`, with $r$ the rank
of $f$ and derived indices $a_1, \ldots, a_r$ and $j$. $\partial_j b$, entry
$j$ of the last axis of the gradient of $b$ on $P$, has the value shape of
$b$ and the free index $j$ of the extent $n$, the dimension of the ambient
coordinate domain of $P$. Empty component-axis lists use the factory
identities: scalar `f[K]` is `f`, and a retained scalar gradient component is
`gradient(b, P)[j]`. $\partial_j$ applies the first row that matches:

| Expression $b$ | $\partial_j b$ |
| -------------- | -------------- |
| $b$ with no domain | `Zero`, with the value shape and the argument spaces of $b$, and the free indices of $b$ followed by $j$ |
| `Coordinate(P)`, $P$ of the dimension of its ambient coordinate domain | `component_tensor(I[i, j], (i,))`, with $I$ the identity matrix of value shape `(n, n)` and a derived index $i$ |
| A terminal, a `Gradient`, a `PushForward`, a node that the identity of Piola matches, or a `PullBack` by a cell or a `SideMap` | `component_tensor(gradient(b, P)[a_1, …, a_s, j], (a_1, …, a_s))`, with $s$ the rank of $b$, so the gradient remains |
| `PullBack(g, φ)`, φ a `DomainMap` | `component_tensor(index_sum(H[A + (l,)] * J[l, j], (l,)), A)`, with `H = pull_back(gradient(g, Q), φ)`, `J = jacobian(φ)`, and $Q$ the ambient domain of φ |
| `Action(A, p, g)` | `Action(∂_j A, p, g)`, with $g$ unchanged |
| Any other expression | The row of [Derivatives of the node classes](#derivatives-of-the-node-classes), with $\delta b = \partial_j b$ |

- In the chain-rule row, `A` contains one derived index per value axis of
  `g`, and `l` is a distinct derived contraction index. The result retains
  the value shape of `g`, with `j` free. For scalar `g`, `A == ()`.
- An operand is constant for the pass if the operand has no domain.
- The gradient of a terminal remains for the library of the terminal, so the
  Hessian of a function is a gradient of a gradient.
- On a mesh $P$ of lower dimension than its ambient coordinate
  domain, `Gradient` is the tangential gradient (see [`Gradient`](expressions.md#gradient)),
  which obeys the same rows. The `Coordinate` row does not match, so
  `gradient(coordinate(P), P)`, the projection $\Pi$, remains, and a
  pull-back by a cell evaluates it.
- A gradient that remains is rewritten on the reference cell by
  [`expand_pull_backs`](#expanding-pull-backs), and the gradients that
  `expand_pull_backs` creates, such as the Jacobian of a cell, are rewritten
  in the next round of the lowering stage.

### Expanding pull-backs

| `expand_pull_backs(f)` |  |
| ---------------------- | - |
| Evaluates              | `PullBack` |
| Input grammar          | $N_3$, in the lowering stage |
| Output grammar         | $N_3 \setminus R$, with $R$ the grammar patterns of the rows below that rewrite |

The pass rewrites every `PullBack(f, ψ)`, inner pull-backs first, with the
first row that matches, and applies the rows again to each `PullBack` that a
rewrite creates. $P$ and $Q$ are the parameter domain and the ambient domain
of $\psi$. For $\psi$ a cell $D_j$ = `cell_map(D, j, e)` of a mesh $D$ with
the topology $T$, $\hat f$ is `pull_back(f, cell_map(T, j, e))` and $m_\psi$
is the law $m$ on $\psi$ (see [`TransformationLaw`](geometry.md#transformationlaw)). For
$\psi$ = `cell_map(Y, i, q)` for a topology or mesh $Y$ and a side $s$,
`B_selected` and
`A_selected` are the maps of [Restriction](geometry.md#restriction), selected by `q`.

| `PullBack(f, ψ)` | Rewrite |
| ---------------- | ------- |
| `Action(A, p, g)` | `Action(PullBack(A, ψ), p, g)`, with $g$ unchanged |
| `Form(e, arguments)` | Make the form fresh for $\psi$, obtaining `Form(e_new, arguments_new)`; return `Form(PullBack(e_new, ψ), arguments_new)` |
| $\psi$ a `SideMap`, $f$ with a domain | Remains until a pull-back by a cell of $P$ |
| $f$ with no domain | $f$, with [type completion](#type-compatibility-of-rewrites) to retain the free indices of $\psi$ |
| `Coordinate(Q)` | The parameterisation of $\psi$ |
| `PullBack(g, SideMap(X, Y, s))`, $\psi$ = `cell_map(Y, i, q)` | `PullBack(PullBack(g, B_selected), A_selected)` |
| `PullBack(g, φ)`, φ a `SideMap` or a cell | Remains, so a pull-back nests only inside a read on a cell |
| `PullBack(g, φ)`, φ a `DomainMap` | `PullBack(g, domain_map(P, pull_back(φ.parameterisation, ψ), φ.ambient_domain))`, the pull-back by the image of $\psi$ under φ (see [`PullBack`](geometry.md#pullback)) |
| `PushForward(f, D, m)`, $\psi$ a cell of $D$ | $m_\psi(\hat f)$ |
| `Divergence(PushForward(f, D, m), D)` of the identity of Piola, $\psi$ a cell of $D$ | $m'_\psi$ of `divergence(f̂, K_j)`, with $m'$ = `transformation_law((), 1)` |
| `Curl(PushForward(f, D, m), D)` of the identity of Piola, $\psi$ a cell of $D$ | $m'_\psi$ of `curl(f̂, K_j)`, with $m'$ = `transformation_law(("contravariant",), 1)` for $n = 3$ and `transformation_law((), 1)` for $n = 2$ |
| `PushForward`, `Divergence` or `Curl` | Remains, read at the points of $P$ |
| `Gradient(g, Q)`, $\psi$ a cell of $Q$ | $m_\psi(\nabla_P(\psi^* g))$, with $m$ = `transformation_law(("covariant",), 0)` |
| `Gradient(g, Q)` | Remains, the gradient on $Q$ read at the points of $P$ |
| `Normal(SideMap(D, E, s))`, $\psi$ = `cell_map(E, i, q)` | $w / \lvert w \rvert$ of [`Normal`](geometry.md#normal) |
| A terminal | Remains, for the library of the terminal |
| Any other node class | The node of the pull-backs of the operands that are expressions |

- The `Form` row pulls back only the body, retaining the declarations after
  capture avoidance (see [Substitution and binding](#substitution-and-binding)).
  Its argument spaces are unchanged. Substitution still visits the body and
  declarations according to their binding scopes.
- The `Action` row transforms only the form's output. The last row moves the
  pull-back into the expression operands of the other unmatched node classes,
  with a binding node made fresh for $\psi$ first.
- A restriction is a pull-back by a `SideMap` (see
  [Restriction](geometry.md#restriction)), so this pass also lowers restrictions. A
  pull-back by a `SideMap` that no pull-back by a cell encloses, such as on an
  expression on $P$ that is not integrated, remains.
- The composition of maps lowers the parameterisation of the map that the
  composition builds as a separate expression root (see
  [Evaluation traversal](#evaluation-traversal)). Composition through a
  `CellMap` remains nested, preserving its explicit cell selection.
- The rows of the identity of Piola read no derivative of $J$. A
  `PushForward` and a node of the identity of Piola under a map that is not a
  cell remain, as a function on a mesh does.
- Only a pull-back by a cell moves a gradient on a mesh to the
  reference cell, so the pass inverts the Jacobian of a cell only. The
  gradient of a pull-back by another map is the chain rule of
  [`apply_gradients`](#applying-gradients).
- Remaining cell reads are evaluated under the library contract in
  [Functions in a space](forms.md#functions-in-a-space).

## The core normal form

The core normal form is $N_4 = N_3 \setminus (G \cup R)$, the output contract of
the standard lowering recipe. Its expression nodes, as visited by the evaluation
traversal, are instances of concrete core classes or downstream terminals.
The following table lists all allowed node classes and their conditions:

| Chapter | Node class | Allowed in the core normal form |
| ------- | ---------- | ------------------------------- |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Integer`, `Float`, `Complex`, `Zero` | Always |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Sum`, `Product`, `Division`, `Power`, `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm` | Always |
| [Tensors and indexing](expressions.md#tensors-and-indexing) | `Tensor`, `Indexed`, `ComponentTensor`, `IndexSum` | Always, including index notation that a rewrite could reduce, such as an `Indexed` of a `ComponentTensor`, which the compiler rewrites (see [Evaluation](#evaluation)) |
| [Tensor operations](expressions.md#tensor-operations) | `Inverse`, `Determinant` | Always, for the compiler to evaluate |
| [Calculus on a domain](expressions.md#calculus-on-a-domain) | `Coordinate` | Always |
| [Calculus on a domain](expressions.md#calculus-on-a-domain) | `Gradient` | Of a terminal, a `Gradient`, a `PushForward`, a node of the identity of Piola, or a `PullBack` by a cell or a `SideMap`, other than `Coordinate(P)` on its own domain $P$ of the dimension of its ambient coordinate domain |
| [Calculus on a domain](expressions.md#calculus-on-a-domain) | `Divergence`, `Curl` | Of the identity of Piola, under no pull-back by a cell |
| [Calculus on a domain](expressions.md#calculus-on-a-domain) | `Integral` | Over a domain that is not a `Mesh`, such as a reference cell or `RD(0)` |
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

Every expression that a pass builds from domain maps, such as the Jacobian,
the inverse Jacobian, the density and the normal of a side, depends only on
those domain maps and remains a subexpression of the result. A form compiler
can therefore evaluate such an expression once for the values that the domain
maps read, such as the coordinates of a mesh, and use the value in every
integral that contains the expression.

In particular, the core normal form may still contain:

- A pull-back by a `SideMap` that no pull-back by a cell encloses, such as a
  restriction that is not integrated.
- A `Gradient` of the coordinates of a mesh of lower dimension
  than its ambient coordinate domain, the projection $\Pi$, and a
  `Gradient`, a `PushForward` or a node of the identity of Piola on a mesh
  not yet read through a cell map, including under an unexpanded side restriction.
