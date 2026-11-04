# Rewriting and lowering

Rewriting replaces expressions without changing their mathematical value.
Lowering applies rewrites to meet the requirements of later computations.

Rule IDs stay with their rules when moved. They cover rewrite rows, including
retained cases, and shared type, traversal and substitution rules. Shared conditions and notation apply to every row. Each rewrite table gives
the input pattern and its conditions, then the replacement or retained result. Contracts and
recipe schedules are not numbered. New rules take new IDs; existing IDs are
not renumbered or reused. The prefixes identify the groups below; numbers follow
reading order within each group at this renumbering. The
[renumbering record](rule-renumbering.md) maps the former IDs to these IDs.

## Shared requirements

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SH01 | Every rewrite | For a valid input expression, a rewrite preserves its mathematical value and [binding](#substitution-and-binding) and must have a [compatible type](#type-compatibility-of-rewrites). It need not preserve floating-point results or simplify algebra. The [mathematical admissibility](rules.md#mathematical-admissibility) convention applies. Zero times an undefined value remains undefined. |

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

### Type compatibility of rewrites

A replacement $r$ must preserve the value of $e$ on its original domain and
satisfy these conditions.

| Rule | Attribute | Compatibility |
| ---- | --------- | ------------- |
| SH02 | Value shape, scalar type, ordered argument spaces | Equal |
| SH03 | Free indices | Same names and extents, up to order and consistent bound-name renaming |
| SH04 | Domain | Equal, or $r$ has no domain and its [constant-field interpretation](expressions.md#domain-free-expressions) equals $e$ on the original domain |

For example, a constant pull-back or the gradient of coordinates may lose its
domain. Code that needs the original domain must keep it separately.

*Type completion* restores scalar type and free indices lost in a rewrite.
Apply it before returning a table's result `r` or rebuilding its parent.
For the replaced expression `e`, let

```python
z = Zero(e.value_shape, e.free_indices, e.argument_spaces, e.scalar_type)
```

Apply both steps in order. SH05 uses tensor arithmetic; SH06 repeats once
for each missing index, in the order of `e.free_indices`.

| Rule | Condition | Completion |
| ---- | --------- | ---------- |
| SH05 | `e` complex, `r` real | `Complex(Integer(1), Integer(0)) * r` |
| SH06 | `r` lacks index `i` of extent `n` from `e.free_indices` | `component_tensor(indexed(broadcast_to(r, (n,) + r.value_shape), (i,) + K), K)`, with `K` one derived index per value axis of `r`, fresh for `e`, `r` and `i` |

For scalar `r`, `K == ()` and the result is `indexed(broadcast_to(r, (n,)), (i,))`.
For tensors, binding `K` restores the value shape. These selections retain the
missing indices without adding arithmetic. `z` is still used for actual zero
results.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SH07 | Result before type completion | The raw result must already have the required shape and argument spaces, matching extents for retained indices and no extra free indices. A real `e` requires a real result. A row giving zero uses `z`; this is what the tables mean by zero with the type of an expression. |
| SH08 | Promotion and completion order | `Action` and `DirectionalDerivative` also promote real values before substitution for complex functions, member by member for tuples; `replace` does not. Promotion must precede evaluation of, for example, a square root at a negative real point. Index completion must precede rebuilding a binder. |

### Evaluation traversal

| Rule | Dependency | Evaluation |
| ---- | ---------- | ---------- |
| SH09 | Root and value operands | Visit recursively in construction order; tuples member by member, including both `Action` operands |
| SH10 | Binding declarations; non-expression `Base` objects | Do not visit |
| SH11 | Domains, maps, spaces, indices and labels | May inspect, but do not traverse or modify stored expressions |
| SH12 | Expressions exposed by a rule | Visit as separate roots; includes newly built map parameterisations |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SH13 | Prepared geometry | Geometry is *prepared* when it satisfies core expressions, actions expanded, derivatives expanded and mesh integrals expanded. Reference recipes require the supplier to prepare geometry before exposing it. Constructors do not require this; another strategy may process unprepared geometry with suitable rules. |

Normal-form checks visit these same nodes. Stored dependencies are checked
separately when needed; labels need no lowering. Substitution follows more
dependencies, as described below.

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
| Tensor arithmetic expanded | Every `Sum`, `Product` and `Division` has scalar value shape |
| Reference integration | Derivatives expanded and mesh integrals expanded |
| Core normal form | The allowed patterns in [The core normal form](#the-core-normal-form) |

A check reports which nodes fail a contract and why, without changing the
expression. Preservation requirements cover domains and operators. Later steps,
including quadrature and tabulation, have their own requirements.

### Notation

| Symbol | Meaning |
| ------ | ------- |
| $\delta b$ | The derivative of $b$; each differentiation rule defines its target |
| $\mathrm{D}b$ | `DirectionalDerivative(b, x, a, h)`; $x$, $a$ and $h$ are its variable, point and direction |
| $\partial_j b$ | Component $j$ of the gradient of $b$ on $P$, with $b$'s value shape and an additional free index $j$ of ambient extent $n$ |
| $P$, $Q$ for `PullBack(f, ψ)` | The parameter domain and ambient domain of $\psi$ |
| Cell | A `CellMap` of a topology or mesh |
| $k$, $n$ in mapping rules | The topological dimension and ambient coordinate dimension of the mapped mesh |
| $A^+$ | `inverse(A)`; rectangular operands satisfy the requirements of [Inverse](expressions.md#inverse) |

| Table entry | Meaning |
| ----------- | ------- |
| Commuting operation | Apply the operation to the listed operands, then rebuild the same node with all other operands unchanged |
| Retained result | Keep the matched expression unchanged for a later rule or the supplying library |
| Excluded pattern | The rule does not apply, even if the rule for that pattern is not selected |
| Derivative identity | Applies where the derivative and the quantities used by the identity are defined, under mathematical admissibility |

### Independence

Independence concerns the operation being applied, not an expression in isolation.
The reference recipes use the checks below.

| Rule | Operation or question | Meaning or check |
| ---- | --------------------- | ---------------- |
| SH14 | Mathematical independence | Varying the operation's target leaves the operand's value unchanged, with other inputs fixed, wherever the expressions are defined |
| SH15 | Directional differentiation with variable $x$ | No member of $x$ occurs free in the operand under the value walk SB22; a single function counts as one member |
| SH16 | Spatial differentiation on $P$ | The operand has domain `None`; interpret it on $P$ by [Domain-free expressions](expressions.md#domain-free-expressions) |
| SH17 | `Action(f, p, g)` | Each unchanged expression operand contributes none of the selected argument slot, according to the node's argument-space definition. Slots are identified by position, not by equality of spaces. Binding declarations and fixed non-expression data do not contribute slots. |
| SH18 | Stronger independence checks | A strategy may use an explicit analysis proving SH14. Reference recipes use SH15–SH17 only. Failure of a check means independence is unknown, not that dependence is proved. Any constant replacement must preserve mathematical admissibility and type. |

| Example | Consequence |
| ------- | ----------- |
| `u - u` | Still contains a free `u`; cancellation requires a separate proof or rewrite |
| `integral(u, P)` | Has no spatial domain, but a free `u` remains for directional differentiation |
| Two argument slots with the same space | Remain distinct slots for an action |

### Commuting rules

Differentiating a sum gives a sum of derivatives. This and the other commuting
rules use SH15–SH17 for unchanged expression operands. Binding declarations
and fixed non-expression data are not differentiated. Specific rules take
precedence. Rename bound names where needed to avoid capturing names in inserted operands (see [Substitution and binding](#substitution-and-binding)).

The commuting rows are grouped under [Structural rules](#structural-rules)
and [Linearity](#linearity). Directional differentiation uses a real parameter
and holds domains fixed. `RealPart` and `ImaginaryPart` have no argument spaces.

## Substitution and binding

Substitution must respect local names. For example, replacing a free `i` must
not change an `i` bound by an index sum.

### Free occurrences

A name is an `Index` or a `Function`. Free occurrence uses the same binding
rules for two walks, chosen by the operation.

| Rule | Walk or condition | Definition |
| ---- | ----------------- | ---------- |
| SB20 | Free occurrence of a name | Encounter a structurally equal name on the chosen walk with no enclosing binder binding it in that operand. Binding declarations do not themselves count as free occurrences. Names mentioned only in labels do not count. |
| SB21 | Construction walk, used by substitution | Follow construction dependencies, including those stored in domains, maps and spaces, respecting the scopes SB02–SB04. Label reconstruction follows [Derived indices](expressions.md#derived-indices), but label contents add no free occurrences. |
| SB22 | Value walk, used by directional differentiation | Follow expression value operands SB01, respecting SB02–SB04. A `Function` contributes itself; do not enter its space or label. Do not enter domains, maps or other stored non-expression dependencies. |

Tuples are walked member by member. A name may be bound in one operand and
free in another. Unless stated otherwise, substitution and capture avoidance
use the construction walk; directional differentiation uses the value walk.

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

To make names fresh for $s_1, \ldots, s_n$, let `body` be the tuple of operands
in the binder's scope. Use these names.

| Rule | Bound name | Fresh name |
| ---- | ---------- | ---------- |
| SB11 | Index `i` | `Index(("uflx.rename", i, body, (s_1, ..., s_n)))` |
| SB12 | Function `v` in $V$ | `Function(V, ("uflx.rename", v, body, (s_1, ..., s_n)))` |

The [derived labels](expressions.md#derived-indices) keep these names distinct
and absent from `body` and the inserted expressions.

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

In these examples, `j` is a fresh index.

| Substitution | Result |
| ------------ | ------ |
| `replace(index_sum(a[i] * c, (i,)), {c: b[i]})` | `index_sum(a[j] * b[i], (j,))`; `i` remains free |
| `replace(form(a, (v,)), {v: w})` | The original form; `v` is bound |
| `replace(DirectionalDerivative(f, u, u, h), {u: u0})` | `DirectionalDerivative(f, u, u0, h)`; the point is outside the binding scope |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SB15 | Compositions | Substitution commutes with a fixed composition when the static-type choice is unchanged. For example, replacing `A` by `B` in `transpose(A)` gives `transpose(B)`, with the derived indices of `B`. |

### Applying forms and preserving binding

An action supplies a value for a form's argument. These rules cover explicit
forms and binding during actions, derivatives and pull-backs.

| Rule | Form $f$ | Rewrite of `Action(f, p, g)` |
| --- | --- | --- |
| SB16 | `Form(e, arguments)`, $p \ge n$ with $n =$ `len(e.argument_spaces)` | Remove entry $p-n$, leaving `others`. Substitute $g$ for that entry in $e$, obtaining `e_new`; for `conjugate(v)`, substitute `conjugate(g)` for `v`, and for tuple entries substitute member by member simultaneously. Return `Form(e_new, others)`, or `e_new` if no entries remain. Preserve binding by SB07–SB14. |

| Rule | Applies to | Requirement |
| --- | --- | --- |
| SB17 | Binding during action substitution | Substitution preserves all enclosing and inserted binders. In particular, commuting through `DirectionalDerivative` makes its variable fresh for `g`; occurrences in `g` are not differentiated. |
| SB18 | Derivative substitution | For `DirectionalDerivative(b, x, a, h)`, match original operands before substituting `[x := a]` into undifferentiated operands; insert `h` unchanged. This substitution holds domains, maps and function spaces fixed, as specified by [DirectionalDerivative](geometry.md#directionalderivative). |

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| SB19 | `Form(e, arguments)` | Make the form fresh for $\psi$, obtaining `Form(e_new, arguments_new)`; return `Form(PullBack(e_new, ψ), arguments_new)`, preserving argument spaces |

## Structural rules

These rules move operations through expression structure or leave them for
another rule or library. Commuting uses the shared conditions above.

| Rule | Node class | Operands | Operations that commute |
| --- | --- | --- | --- |
| ST01 | `Indexed`, `ComponentTensor`, `IndexSum`, `BroadcastTo`, `Form`, `DirectionalDerivative`; for `Action` on a `Form`, require $p<n$, where $n =$ `len(e.argument_spaces)`; nested derivatives follow D13 | The expression operand (`value` for `BroadcastTo`, `f` for `DirectionalDerivative`, `e` for `Form`) | `Action`, `DirectionalDerivative`, spatial partial derivative |

| Rule | Form $f$ | Rewrite of `Action(f, p, g)` |
| --- | --- | --- |
| ST02 | `Product(a, b)`, with $n$ the arity of $a$ | `Product(Action(a, p, g), b)` for $p < n$, otherwise `Product(a, Action(b, p - n, g))` |
| ST03 | `Gradient`, `PullBack` | Remains; the operation transforms the form's output |
| ST04 | `Function` or a downstream terminal | Remains, for the supplying library |
| ST05 | `Action` | Remains, as a further application of a remaining action |
| ST06 | A core expression with `argument_spaces == ()` | No valid `Action` can select an argument position |

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| ST07 | `Gradient(Action(B, p, g), P)` | `Action(Gradient(B, P), p, g)` |
| ST12 | `Gradient(b, P)`, $b$ with a domain and either a terminal, `Gradient`, `PushForward`, node matching CC06 or CC07, or `PullBack` by a cell or `SideMap`; exclude `Coordinate(P)` when $P$ has its ambient dimension | Retain the whole gradient unchanged |
| ST08 | `Gradient(f, P)`, $f$ not an `Action` and the gradient does not match ST12 | `component_tensor(∂_j f[A], A + (j,))`, using the component rules selected by [Applying gradients](#applying-gradients) |
| ST09 | $\partial_j b$, with `Gradient(b, P)` matching ST12 | `component_tensor(gradient(b, P)[B + (j,)], B)`, with one derived index in `B` per value axis of $b$; ST12 retains the introduced gradient |
| ST10 | $\partial_j$ `Action(B, p, g)` | `Action(∂_j B, p, g)`, with $g$ unchanged |

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| ST11 | `Action(A, p, g)` | `Action(PullBack(A, ψ), p, g)`, with $g$ unchanged |

Gradient reconstruction uses derived indices $A$ and $j$, with one index in
$A$ per value axis of $f$. Empty lists use the factory identities. The supplying
library evaluates remaining terminal gradients. Downstream terminals may declare
argument spaces, for example for assembled vectors or matrices; their libraries
evaluate the remaining actions.

### Tensor arithmetic components

These rules expand tensor arithmetic into scalar operations. Each applies only
to a non-scalar result. `K` has one derived index per result axis, using
`Index(("uflx.expand_arithmetic", f, m))` for node `f` and axis `m`.
Choose fresh names under SB07–SB14; operand free indices remain free.

| Rule | Node | Rewrite |
| ---- | ---- | ------- |
| ST13 | `Sum(a, b)` | `component_tensor(Sum(a[K], b[K]), K)` |
| ST14 | `Product(a, b)`, scalar `a` | `component_tensor(Product(a, b[K]), K)` |
| ST15 | `Product(a, b)`, scalar `b` | `component_tensor(Product(a[K], b), K)` |
| ST16 | `Division(a, b)` | `component_tensor(Division(a[K], b), K)` |

## Linearity

These rows apply only to the listed operations and operands, under the shared
commuting conditions. The `Tensor` case constructs a tensor from the results
for its components. Conjugation uses a real differentiation parameter.

| Rule | Node class | Operands | Operations that commute |
| --- | --- | --- | --- |
| L01 | `Sum`, `Tensor` | Every operand | `Action`, `DirectionalDerivative`, spatial partial derivative |
| L02 | `Division` | The numerator | `Action`, `DirectionalDerivative`, spatial partial derivative |
| L03 | `Gradient`, `Divergence`, `Curl`, `PushForward`, `PullBack`, `Integral` | The expression | `DirectionalDerivative`; also `Action` except on `Gradient` or `PullBack` (ST03) |
| L04 | `Conjugate`, `RealPart`, `ImaginaryPart` | The expression | `DirectionalDerivative`, spatial partial derivative |

| Rule | Form $f$ | Rewrite of `Action(f, p, g)` |
| --- | --- | --- |
| L05 | `Conjugate(e)` | `conjugate(action(e, p, conjugate(g)))`, with conjugation member by member for a `tuple` $g$ |
| L07 | `Zero` | Zero with the type of the action |

| Rule | Expression $f$ | Rewrite of the derivative at $a$ |
| --- | --- | --- |
| L10 | $f$ passes the directional independence check SH15 for $x$ | Zero with the type of the derivative |

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| L11 | $\partial_j b$, $b$ passes the spatial independence check SH16 | Zero with $b$'s shape and argument spaces, and its free indices followed by $j$ |

## Differentiation

### Derivatives of the node classes

Directional and spatial derivatives share these rules. Their specific rules
take precedence. The linear derivative rules are in [Linearity](#linearity).
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
| D08 | `Determinant(A)`, $A$ invertible | $\det(A)\sum_{i,j}(A^{-1})_{ij}\,\delta A_{ji}$, with `inverse`, `index_sum` and derived indices $i,j$ |
| D09 | `Inverse(A)`, equal axis extents, including symbolic extents, and $A$ invertible | `-inverse(A) @ delta_A @ inverse(A)`, with `delta_A` = $\delta A$ |
| D10 | `Inverse(A)`, $A$ real with shape `(n, k)`, $k<n$, and full column rank | $-A^+ \, \delta A \, A^+ + A^+ (A^+)^T (\delta A^T - \delta A^T A A^+)$, with `@`, `transpose` and `-` |

Literal values in D19 and D20 are compared exactly, including signed zero;
`Complex` has that value when its imaginary part is zero. These cases apply
before substituting the differentiation point. A variable exponent that happens
to be zero or one at that point still uses D04. An undefined base remains
undefined under mathematical admissibility.

The constant-exponent power rule avoids a logarithm, permitting negative real
bases where differentiable. Results are core expressions.

### Directional derivatives

| Rule | Expression $f$ | Rewrite of the derivative at $a$ |
| --- | --- | --- |
| D11 | $f=x$ for a single function, or $f=x_k$ for a tuple variable | $h$ for a single function, or $h_k$ for the matching tuple member |
| D12 | `Action(b, p, g)` | `Action(Db, p, g) + Action(b, p, Dg)`, with undifferentiated operands evaluated at $a$ |
| D13 | `DirectionalDerivative` | Defer until the inner derivative is expanded |

| Rule | Applies to | Requirement |
| --- | --- | --- |
| D15 | Dependence on the derivative variable | Use the directional independence check SH15 and value walk SB22. Tuple variables vary member by member, including occurrences under conjugation. Dual directions and substitution points may be symbolic forms whose actions the recipe evaluates. |

For real scalar functions `u` and `a` in a space of constants,
`directional_derivative(u * u, u, a, u)` rewrites to `u * a + a * u`.

For a real scalar function `u` in a space of constants and a constant array
`a` of shape `(m,)`, `directional_derivative(u, u, a[i], Integer(1))` lowers
to `indexed(broadcast_to(Integer(1), (m,)), (i,))`. Its value is one for every `i`.
Type completion retains indices even from unused members of tuple points
and directions.

### Spatial derivatives

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| D16 | `Divergence` or `Curl` not matching CC06 or CC07 | Its composition of gradients from [Divergence and Curl](geometry.md#divergence-and-curl) |
| D17 | $\partial_j$ `Coordinate(P)`, $P$ of its ambient dimension | `component_tensor(I[i, j], (i,))`, with $I$ the identity of shape `(n, n)` and a derived $i$ |

On embedded meshes these rules apply to the tangential gradient. The coordinate
gradient remains a projection until the pull-back to a cell.

## Pull-back algebra

Pull-backs express values in a map's parameter coordinates. These rules handle
pointwise operations, coordinates and map composition. Binding is handled by
the substitution rules; restrictions and transformation laws have their own groups.

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| PB01 | $f$ with no domain, neither an `Action` nor a `Form` | The same spatially constant value $f$ under [Domain-free expressions](expressions.md#domain-free-expressions), with [type completion](#type-compatibility-of-rewrites) to retain the free indices of $\psi$ |
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
the [geometry preparation contract](#evaluation-traversal).

## Change of coordinates

### Gradients and mapped values

For a mesh cell $\psi = D_j =$ `cell_map(D, j, e)` with topology $T$, let
$K_j$ be its reference cell, $\hat f$ = `pull_back(f, cell_map(T, j, e))`,
and $m_\psi$ be the law $m$ on $\psi$ (see
[`TransformationLaw`](geometry.md#transformationlaw)).

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| CC01 | $\partial_j$ `PullBack(g, φ)`, φ a `DomainMap` with parameter domain $P$ | `component_tensor(index_sum(H[A + (l,)] * J[l, j], (l,)), A)`, with `H = pull_back(gradient(g, Q), φ)`, `J = jacobian(φ)`, and $Q$ the ambient domain of φ |

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| CC02 | `PushForward(f, D, m)`, $\psi$ a cell of $D$ | $m_\psi(\hat f)$ |
| CC03 | `PushForward`, `Divergence` or `Curl` not matching CC02, CC08 or CC09; exclude R01 | Remains, read at the points of $P$ |
| CC04 | `Gradient(g, Q)`, $\psi$ a cell of $Q$ | $m_\psi(\nabla_P(\psi^* g))$, with $m$ = `transformation_law(("covariant",), 0)` |
| CC05 | `Gradient(g, Q)` not matching CC04; exclude R01 | Remains, the gradient on $Q$ read at the points of $P$ |

In the gradient chain rule, `A` indexes the value axes of `g` and `l` is a
distinct contraction index; scalar `g` uses `A == ()`. Only cell transport
inverts a Jacobian; other maps use the chain rule. Libraries evaluate remaining
cell reads under [Functions in a space](forms.md#functions-in-a-space).

### Piola identities

The retained divergence and curl patterns take precedence over generic
expansion. The cell pull-back rules use no derivative of $J$.

| Rule | Expression or operation | Rewrite |
| --- | --- | --- |
| CC06 | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)` | Remains for the Piola rule of a cell pull-back |
| CC07 | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$ | Remains for the Piola rule of a cell pull-back |

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| CC08 | `Divergence(PushForward(f, D, m), D)`, `m = transformation_law(("contravariant",), 1)`, $\psi$ a cell of $D$ | $m'_\psi$ of `divergence(f̂, K_j)`, with $m'$ = `transformation_law((), 1)` |
| CC09 | `Curl(PushForward(f, D, m), D)`, `m = transformation_law(("covariant",), 0)`, $k=n$, $\psi$ a cell of $D$ | $m'_\psi$ of `curl(f̂, K_j)`, with $m'$ = `transformation_law(("contravariant",), 1)` for $n = 3$ and `transformation_law((), 1)` for $n = 2$ |

## Restriction

For $\psi =$ `cell_map(Y, i, q)` of a topology or mesh $Y$ and a side $s$,
`B_selected` and `A_selected` are the maps of
[Restriction](geometry.md#restriction), selected by `q`.
Cell maps retain explicit selection; unselected side restrictions remain.

| Rule | Operand of `PullBack(·, ψ)` and conditions | Rewrite |
| --- | --- | --- |
| R01 | $\psi$ a `SideMap`, $f$ with a domain and neither an `Action` nor a `Form` | Remains until a pull-back by a cell of $P$ |
| R02 | `PullBack(g, SideMap(X, Y, s))`, $\psi$ = `cell_map(Y, i, q)` | `PullBack(PullBack(g, B_selected), A_selected)` |
| R03 | `Normal(SideMap(D, E, s))`, $\psi$ = `cell_map(E, i, q)`; the cell Jacobian has full column rank and $w \ne 0$ | $w / \lvert w \rvert$, with $w$ defined in [The pull-back of Normal](geometry.md#the-pull-back-of-normal) |

## Integration

Mesh integration becomes a sum over cells, with each cell integral expressed
on its reference cell. The density accounts for the change of coordinates.

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
| I01 | `Integral(f, D)`, $D$ a `Mesh` | `component_tensor(sum_j, K)`, where `sum_j` is the sum over groups $j$ of `index_sum(integral(g_j * density(D_j), K_j), (e_j,))` |

$g_j$ is `pull_back(f_K, piece_j)` when $f$ has domain `D.ambient_domain`,
and `pull_back(f_K, D_j)` otherwise, including domain-free $f$.
The pull-back introduces the internal cell index; `IndexSum` sums over it and
retains every free index of $f$ and the component indices `K`. The outer
`component_tensor` restores the value shape, and is the identity for scalar
$f$. The remaining integrals have coordinate domains.

Expanding the same input again gives the same graph. Changing the order of
substitution and expansion may change bound names when renaming is needed.
The value and free indices must stay the same.

## Pass recipes

All recipes use the shared requirements. Selection rows are tried from top to bottom, and rule IDs within a row from
left to right, including retained cases. Binding and dependence rows state
requirements used throughout the recipe. Fallback selections apply
only when the preceding rules do not match. Rules keep their conditions when
selected by another strategy.

### Visiting and revisiting

These rules apply to the reference recipes. Retention is a rule result, not a
permanent annotation on a node.

| Rule | Event | Required action |
| ---- | ----- | --------------- |
| SH19 | Visit a node | Process its value operands in construction order under SH09–SH11, then rebuild if they changed and select a rule for the resulting node |
| SH20 | A rule replaces a node | Apply type completion, then process the replacement, including its root and any new or changed subexpressions that this recipe handles. Reuse already processed, structurally unchanged subexpressions. Substitution itself still follows SB13. |
| SH21 | A rule retains a node | Finish this visit without reconstructing it through a fallback or scheduling it again unchanged. Its operands have already been processed. Reconsider it if an operand changes or another recipe visits it. |
| SH22 | A rule exposes work for another recipe | Leave that work for the stated recipe. After it changes an operand, rebuild and reconsider affected ancestors. Exposed geometry follows SH12–SH13. |
| SH23 | A recipe stops | No work belonging to that recipe remains. Check its output contract. An unchanged replacement is not progress; an unmatched required operation, repeated state or resource limit is not successful completion. |
| SH24 | A replacement contains $\partial_j b$ | Resolve this component-derivative request with the component rows of `apply_gradients`. It is notation, not a core node. ST08 creates these requests; ST09 ends a request at a retained gradient. No unresolved request may appear in the returned expression. |

A retained node does not hide its operands. For example, a `PullBack` retained
around a gradient must be reconsidered if a later gradient pass changes that
operand into arithmetic. A recipe may revisit an unchanged retained node in a
later round, but its retention rule returns it unchanged without new indices.

### Applying actions

| `apply_actions(f)` |  |
| ------------------ | - |
| Evaluates          | `Action` |
| Input contract     | Core expressions |
| Output contract    | Actions expanded |
| Strategy | SH19–SH23; inner actions first; first matching rule; process new actions, including those exposed by substitution |

| Selection | Rules |
| --------- | ----- |
| Apply explicit forms | SB16 |
| Conjugation and products | L05, ST02 |
| Keep spatially transformed forms | ST03 |
| Commute into operands | L01, ST01, L02–L03 |
| Zero and remaining applications | L07, ST04–ST06 |
| Binding | SB01–SB15, SB17, SB20–SB21 |

### Applying derivatives

| `apply_derivatives(f)` |  |
| ---------------------- | - |
| Evaluates              | `DirectionalDerivative` and actions created by its evaluation |
| Input contract         | Core expressions; actions expanded |
| Output contract        | Actions and derivatives expanded |
| Strategy | SH19–SH23; inner derivatives first; first matching rule; process new derivatives and run `apply_actions` on exposed actions, then reconsider changed ancestors |

| Selection | Rules |
| --------- | ----- |
| Substitution and dependence | SB18, D15, SH15, SB20, SB22 |
| Independent expression | L10 |
| Variable and action derivatives | D11–D12 |
| Defer an inner derivative | D13 |
| Derivative fallback | D01–D02, D19–D20, D03–D10, L04, L01, ST01, L02–L03, in order, with $\delta b = \mathrm{D} b$ |

### Expanding integrals

| `expand_integrals(f)` |  |
| --------------------- | - |
| Evaluates             | `Integral` over a mesh |
| Input contract        | Core expressions |
| Output contract       | Mesh integrals expanded |
| Strategy | SH19–SH23; inner mesh integrals first; process mesh integrals exposed by expansion |

| Selection | Rules |
| --------- | ----- |
| Mesh integration | I01 with the notation and index construction in [Integration](#integration) |

### Applying gradients

| `apply_gradients(f)` |  |
| -------------------- | - |
| Evaluates            | `Gradient`, `Divergence` and `Curl` |
| Input contract       | Core expressions; actions and derivatives expanded |
| Output contract      | Actions, derivatives and gradients expanded |
| Strategy | SH19–SH24; inner expressions first; first matching rule; resolve component derivatives and run `apply_actions` on exposed actions, then reconsider changed ancestors |

| Selection | Rules |
| --------- | ----- |
| Move an action outside the gradient | ST07 |
| Keep Piola patterns and terminal or mapped gradients | CC06–CC07, ST12 |
| Expand divergence and curl | D16 |
| Reconstruct the gradient from components | ST08 |
| Constants and coordinates | L11, D17 |
| Keep terminal and mapped gradients | ST09 |
| Chain rule and actions | CC01, ST10 |
| Derivative fallback | D01–D02, D19–D20, D03–D10, L04, L01, ST01, L02–L03, in order, with $\delta b = \partial_j b$ |

Gradients exposed by pull-backs are processed in the next call to `apply_gradients`.

### Expanding pull-backs

| `expand_pull_backs(f)` |  |
| ---------------------- | - |
| Evaluates              | `PullBack` |
| Input contract         | Core expressions; actions and derivatives expanded; prepared geometry when exposed |
| Output contract        | Actions, derivatives and pull-backs expanded |
| Strategy | SH19–SH23; inner pull-backs first; first matching rule; process new pull-backs and run `apply_actions` on exposed actions. Leave new gradients to the next spatial round. |

| Selection | Rules |
| --------- | ----- |
| Actions and binding | ST11, SB19 |
| Keep a side restriction | R01 |
| Constants and coordinates | PB01–PB02 |
| Select the cell on a side | R02 |
| Nested maps | PB03–PB04 |
| Mapped values | CC02 |
| Piola identities | CC08–CC09 |
| Remaining mapped derivatives | CC03–CC05 |
| Side normal | R03 |
| Terminals and pointwise operations | PB05–PB06 |

### Expanding tensor arithmetic

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `expand_arithmetic(f)` | | `f: AbstractExpression` | An equivalent expression with tensor arithmetic expanded |

| `expand_arithmetic(f)` | |
| ---------------------- | - |
| Evaluates | Tensor-valued `Sum`, `Product` and `Division` |
| Input contract | Core expressions |
| Output contract | Tensor arithmetic expanded |
| Strategy | SH19–SH23; operands first; apply ST13–ST16; retain scalar arithmetic |

## The core normal form

The core normal form combines
the core expressions, actions expanded, derivatives expanded, mesh integrals
expanded, gradients expanded, pull-backs expanded and tensor arithmetic expanded
contracts. The table lists the nodes allowed after `lower`.

| Chapter | Node class | Allowed in the core normal form |
| ------- | ---------- | ------------------------------- |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Integer`, `Float`, `Complex`, `Zero` | Always |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Sum`, `Product`, `Division` | Scalar value shape |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Power`, `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm` | Always |
| [Tensors and indexing](expressions.md#tensors-and-indexing) | `Tensor`, `BroadcastTo`, `Indexed`, `ComponentTensor`, `IndexSum` | Always; index simplification is left to the compiler |
| [Tensor operations](expressions.md#tensor-operations) | `Inverse`, `Determinant` | Always, for the compiler to evaluate |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Coordinate` | Always |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Gradient` | Matches ST12 |
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
| Stable output violates core normal form | `ValueError`, identifying offending nodes |
| A complete spatial round returns to an earlier state | `ValueError`; a cycle is not completion |

It runs these passes in order.

| Recipe step | Required result |
| ----------- | --------------- |
| [`apply_actions`](#applying-actions) | Actions expanded |
| [`apply_derivatives`](#applying-derivatives), including action cleanup | Actions and derivatives expanded |
| [`expand_integrals`](#expanding-integrals) | Also mesh integrals expanded |
| Repeat [`apply_gradients`](#applying-gradients), then [`expand_pull_backs`](#expanding-pull-backs), with action cleanup and reconsideration under SH19–SH24 | Stop under the conditions below |
| [`expand_arithmetic`](#expanding-tensor-arithmetic) | Core normal form |

Compare states after complete spatial rounds using structural equality,
including bound names and labels. Fresh naming alone is not progress.

| Round result | Next step |
| ------------ | --------- |
| Neither spatial pass nor action cleanup changed the expression; all contracts except tensor arithmetic expanded hold | Proceed to `expand_arithmetic` |
| Neither changed it; another core normal-form contract fails | Raise `ValueError` with the remaining offending nodes |
| Some step changed it, but the final state equals an earlier round boundary | Report a cycle; do not treat it as an unchanged round |
| Otherwise | Run the next round |

### Progress of the reference recipe

Expression size may grow, so node count is not a progress measure. The rules
remove an operation, move it to operands, or stop at an explicit retained case.

| Rule family | Progress |
| ----------- | -------- |
| Actions and directional derivatives | Apply or differentiate the current node, then process operations on its operands. An inner derivative is processed before an outer one. Retained actions are not expanded again. |
| Mesh integration | Replace a mesh integral by reference-cell integrals. Repeating the pass cannot expand those as mesh integrals. |
| Gradient reconstruction | Apply ST08 only to a gradient outside ST12. Component rules descend through operands. ST09 introduces a gradient that matches ST12, so reconstruction cannot wrap it again. |
| Action transport | ST07, ST10 and ST11 move actions outside spatial operations. ST03 prevents action expansion from reversing this move. |
| Pull-back algebra | Move a pull-back into operands, compose an existing map, or retain a read. Composition exposes the existing map's parameterisation for processing; it does not recreate the original map composition. |
| Cell transport and Piola identities | Replace the matched cell pull-back by its transformation formula. New gradients are handled in the next spatial round and stop at ST12 where appropriate. |
| Tensor arithmetic | ST13–ST16 replace a tensor arithmetic node by scalar arithmetic and index notation. They introduce no tensor arithmetic nodes. |
| Retained cases | Do not create fresh indices or schedule another rewrite of the unchanged node |

The finite, acyclic input and stored geometry supply the operands and map
parameterisations in these steps. New recipes or extension rules must justify
their own progress; the stopping checks above do not prove termination for an
arbitrary rule selection.

### Stability examples

Here `u` is a scalar function on `P` and `v` a vector function on `P`.

| Input | Required behaviour |
| ----- | ------------------ |
| `gradient(u, P)` | ST12 returns it unchanged, without a `ComponentTensor` wrapper |
| `gradient(gradient(u, P), P)` | Visit the inner gradient, then retain both under ST12 |
| `gradient(v, P)` | Retain the whole tensor gradient; do not expand and rebind its components |
| `gradient(u * u, P)` | Expand the product derivative; retain each resulting `gradient(u, P)` |
| Gradient of a domain-free value | ST12 does not apply; L11 gives the typed zero |
| `gradient(coordinate(P), P)` on a full-dimensional domain | ST12 does not apply; D17 gives the identity |
| A retained pull-back whose gradient operand later becomes arithmetic | Reconsider the rebuilt pull-back and apply its arithmetic rules |
| A completed `lower(f)` passed to `lower` again | Return a structurally equal expression |

### Worked rule combinations

These checks compare values, shapes and free indices. They do not require
arithmetic simplification. Zero terms are omitted from displayed results,
and fresh bound names are shown schematically.

#### Nested derivatives

Let `u`, `a`, `b`, `h` and `k` be real scalar functions in a constant space,
with `a`, `b`, `h` and `k` independent of `u`.

| Input | Rule trace and result |
| ----- | --------------------- |
| `directional_derivative(directional_derivative(u * u, u, u, h), u, a, k)` | The inner derivative gives `h * u + u * h` by D01, D11 and SB18. Its point `u` is free under SB04. The outer derivative gives `h * k + k * h`. |
| `directional_derivative(directional_derivative(u * u, u, b, h), u, a, k)` | The inner derivative gives `h * b + b * h`. SH15 then finds no free `u`, so L10 gives the typed zero. |

#### Capture avoidance

Let `u` and `u0` be real scalar functions in a constant space, and `a` and `b`
be constant arrays. The index `i` has the same extent in both arrays.

| Input | Rule trace and result |
| ----- | --------------------- |
| `directional_derivative(index_sum(u * a[i], (i,)), u, u0, b[i])` | SB06–SB09 rename the sum's bound index to fresh `j` before ST01 moves differentiation inside. D01 and D11 give `index_sum(b[i] * a[j], (j,))`. The direction's `i` stays free. |

For `a = (2, 3)` and `b = (7, 11)`, the result has values `(35, 55)` over
free `i`. Capturing `i` would instead give the scalar `47`.

#### Several missing indices

Let `u` be a real scalar function in a constant space, `a` and `b` constant
arrays of extents `m` and `n`, and `i` and `j` distinct indices.

| Step | Result |
| ---- | ------ |
| Differentiate `directional_derivative(u, u, a[i] + b[j], Integer(1))` | D11 gives the raw value `Integer(1)`, missing both indices |
| SH06 for `i` | `r = indexed(broadcast_to(Integer(1), (m,)), (i,))` |
| SH06 for `j` | `indexed(broadcast_to(r, (n,)), (j,))` |
| Bind and sum both indices | The value is `m * n`, including zero when either extent is zero |

For a vector result of shape `(d,)`, SH06 instead uses a fresh component
index `k` and returns
`component_tensor(indexed(broadcast_to(r, (m, d)), (i, k)), (k,))`.
The result keeps shape `(d,)` and adds free `i`; a second missing index is
handled in the same way without capturing `k`.

#### Empty broadcasts

| Input | Required result |
| ----- | --------------- |
| `index_sum(indexed(broadcast_to(square_root(Integer(-1)), (0,)), (i,)), (i,))` | Value zero with shape `()` and no free indices. Neither the empty broadcast nor the empty sum evaluates the real square root. The prototype remains subject to static construction checks. |

This does not make `square_root(Integer(-1))` a defined real value. Symbolic
traversal under SH19 is distinct from evaluating an array member.

#### Piola identity exposed by an action

Let `v` and `u` be vector functions on the topology of mesh `D`, with geometry
independent of them. Set `m = transformation_law(("contravariant",), 1)`,
`F = form(push_forward(v, D, m), (v,))`, and let `ψ` be a cell map of `D`.

| Step | Result |
| ---- | ------ |
| Start with `pull_back(divergence(action(F, 0, u), D), ψ)` | The Piola pattern is hidden by the action |
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

Other recipes may choose individual rules or rule families, provided their
conditions hold and geometry is prepared as required. For example,
[shape differentiation](extensions.md#shape-derivatives) lowers spatially before
and after differentiation. Index-notation simplification, such as reducing an
`Indexed` of a `ComponentTensor`, remains a compiler operation.
