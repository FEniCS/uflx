# Rewrite machinery

## Shared requirements

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SH01 | Every rewrite | For a valid input expression, a rewrite preserves its mathematical value and [binding](substitution.md#substitution-and-binding) and must have a [compatible type](#type-compatibility-of-rewrites). Floating-point evaluation of the original and rewritten expressions may differ. Algebraic simplification is not required. The [mathematical admissibility](rules.md#mathematical-admissibility) convention applies. Zero times an undefined value remains undefined. |

The lazy nodes are `Apply`, `DirectionalDerivative`, mesh `Integral`,
`PullBack`, `Gradient`, `Divergence`, `Curl` and `PushForward`. They store
operations to evaluate later.

### Transformations, recipes and contracts

| Component | Responsibility |
| --------- | -------------- |
| Rule | Match part of an expression, check conditions and replace it; use the shared rules for binding, type completion and rebuilding nodes |
| Recipe | Select and schedule rules, with input, output and preservation requirements |
| Contract | State a requirement on an expression, checked without rewriting it |

The passes in [Differentiation](differentiation.md) and
[Applications and mappings](lowering.md) are reference recipes. Rules cannot
invoke whole passes. No rule objects or global registry are required.

A recipe’s *strategy* is its choice of rules, traversal and stopping conditions.
Contracts can be combined and apply to the nodes visited by the
[evaluation traversal](#evaluation-traversal).

| Condition | Requirement |
| ------------------ | ----------- |
| Specific identity or a pattern that must remain | Takes precedence over generic expansion, even when its rule is not selected |
| Unmatched node | May remain during partial rewriting; required lowering must satisfy its output contract |
| Optional simplification | Need not eliminate any node class |
| Repetition | Must make progress without cycling; reaching a time or step limit does not count as success |
| Expression no longer changes | Does not mean all equivalent expressions have the same structure, or that different recipes give the same output |

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

### Evaluation traversal

| Rule | Dependency | Evaluation |
| ---- | ---------- | ---------- |
| SH09 | Root and value operands | Visit recursively in construction order; tuples member by member, including the form and all supplied arguments of `Apply` |
| SH10 | Binding declarations; non-expression `Base` objects | Do not visit |
| SH11 | Domains, maps, spaces, indices and labels | May inspect, but do not traverse or modify stored expressions |
| SH12 | Expressions exposed by a rule | Visit as separate roots; includes newly built map parameterisations |

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SH13 | Prepared geometry | Geometry is *prepared* when it satisfies core expressions, applications expanded, derivatives expanded and mesh integrals expanded. Reference recipes require the supplier to prepare geometry before exposing it. Constructors do not require this; another strategy may process unprepared geometry with suitable rules. |

Normal-form checks visit these same nodes. Stored dependencies are checked
separately when needed; labels need no lowering. Substitution follows more
dependencies, as described in [Substitution and binding](substitution.md#substitution-and-binding).

### Visiting and revisiting

Reference recipes use the procedure below. Try selection rows top to bottom
and rules left to right; the first match wins, including retained cases.
Binding and dependence requirements apply throughout.

| Rule | Event | Required application |
| ---- | ----- | --------------- |
| SH19 | Visit a node | Process its value operands in construction order under SH09–SH11, then rebuild if they changed and select a rule for the resulting node |
| SH20 | A rule replaces a node | Apply type completion, then process the replacement, including its root and any new or changed subexpressions that this recipe handles. Reuse already processed, structurally unchanged subexpressions. Substitution itself still follows SB13. |
| SH21 | A rule retains a node | Finish this visit without reconstructing it through a fallback or scheduling it again unchanged. Its operands have already been processed. Reconsider it if an operand changes or another recipe visits it. |
| SH22 | A rule exposes work for another recipe | Leave that work for the stated recipe. After it changes an operand, rebuild and reconsider affected ancestors. Exposed geometry follows SH12–SH13. |
| SH23 | A recipe stops | No work belonging to that recipe remains. Check its output contract. An unchanged replacement is not progress; an unmatched required operation, repeated state or resource limit is not successful completion. |

For example, reconsider a retained `PullBack` when gradient expansion changes
its operand into arithmetic. Retaining an unchanged node creates no new indices.

## Shared mathematical rules

### Independence

Independence concerns the operation being applied, not an expression in isolation.
The reference recipes use the checks below.

| Rule | Operation or question | Meaning or check |
| ---- | --------------------- | ---------------- |
| SH14 | Mathematical independence | Varying the operation's target leaves the operand's value unchanged, with other inputs fixed, wherever the expressions are defined |
| SH15 | Directional differentiation with variable $x$ | After geometry preparation D21–D23, no member of $x$ occurs free in the operand under the value walk SB22; a single function counts as one member. This check alone must not discard a construction dependency. |
| SH16 | Spatial differentiation on $P$ | The operand has domain `None`; interpret it on $P$ by [Domain-free expressions](expressions.md#domain-free-expressions) |
| SH17 | `Apply(f, G)` | Each unchanged expression operand has no argument slots, according to the node's argument-space definition. Ordered slots determine how `G` is distributed; equality of spaces never identifies slots. Binding declarations and fixed non-expression data do not contribute slots. |
| SH18 | Stronger independence checks | A strategy may use an explicit analysis proving SH14. Reference recipes use SH15–SH17 only. Failure of a check means independence is unknown, not that dependence is proved. Any constant replacement must preserve mathematical admissibility and type. |

| Example | Consequence |
| ------- | ----------- |
| `u - u` | Still contains a free `u`; cancellation requires a separate proof or rewrite |
| `integral(u, P)` | Has no spatial domain, but a free `u` remains for directional differentiation |
| Two argument slots with the same space | Remain distinct slots for an application |

### Commuting rules

Commuting uses SH15–SH17 for unchanged operands and capture-avoiding renaming
([Substitution](substitution.md)). Specific rules take precedence. Before holding
binding declarations or non-expression data fixed, check derivative dependencies
under D21–D23. Directional differentiation uses a real parameter;
`RealPart` and `ImaginaryPart` have no argument spaces.

### Structural rules

| Rule | Node class | Operands | Operations that commute |
| --- | --- | --- | --- |
| ST01 | `Indexed`, `ComponentTensor`, `IndexSum`, `BroadcastTo`, `Form`, `DirectionalDerivative`; full application of a `Form` uses SB16 instead; nested derivatives follow D13 | The expression operand (`value` for `BroadcastTo`, `f` for `DirectionalDerivative`, `e` for `Form`) | `Apply`, `DirectionalDerivative`, spatial partial derivative |

### Linearity

The commuting conditions apply. `Tensor` collects the component results;
conjugation uses a real differentiation parameter.

| Rule | Node class | Operands | Operations that commute |
| --- | --- | --- | --- |
| L01 | `Sum`, `Tensor` | Every operand | `Apply`, `DirectionalDerivative`, spatial partial derivative |
| L02 | `Division` | The numerator | `Apply`, `DirectionalDerivative`, spatial partial derivative |
| L03 | `Gradient`, `Divergence`, `Curl`, `PushForward`, `PullBack`, `Integral` | The expression | `DirectionalDerivative`; also `Apply` except on `Gradient` or `PullBack` (ST03) |
| L04 | `Conjugate`, `RealPart`, `ImaginaryPart` | The expression | `DirectionalDerivative`, spatial partial derivative |

## Rewrite mechanisms

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

Differentiating a real scalar constant-space function `u` at `a[i]` in direction
one gives one for every `i`. Type completion retains that family index.

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

For scalar `r`, `K == ()`, giving `indexed(broadcast_to(r, (n,)), (i,))`.
For tensors, binding `K` restores the shape. Completion adds no arithmetic;
actual zero results use `z`.

| Rule | Applies to | Requirement |
| ---- | ---------- | ----------- |
| SH07 | Result before type completion | The raw result must already have the required shape and argument spaces, matching extents for retained indices and no extra free indices. A real `e` requires a real result. A row giving zero uses `z`; this is what the tables mean by zero with the type of an expression. |
| SH08 | Promotion and completion order | `Apply`, `partial` and `DirectionalDerivative` promote real values before substitution for complex functions, member by member for tuples; `replace` does not. Promotion precedes evaluation, for example of a square root at a negative real point. Index completion precedes rebuilding a binder. |

### Grammars

A grammar is a set of allowed expressions. A contract check tests membership
of that set (see [New questions](extensions.md#new-questions)).

A grammar pattern inspects a node and objects within a bounded number of links,
without running a pass. It may classify a domain or map without adding that
dependency to the traversal. For a grammar $N$ and patterns $P$, $N \setminus P$
contains the expressions of $N$ with no visited node matching $P$. Combining
contracts means taking the intersection of their grammars.

A recipe maps its input grammar to its output grammar under stated conditions.
A rule guarantees only its replacement and preservation requirements.
Extensions must supply grammars and rules before lowering their nodes to core.

The form layer contains every `Apply` and every node with nonempty argument
spaces. The integrand layer contains all other nodes. A contract may restrict
these separately. A `Function` can satisfy core normal form even if its stored
parameterisation contains a reducible gradient. That parameterisation is checked
separately when exposed.

## Type-completion examples

Displayed results omit zero terms and use schematic fresh names;
rewriting need not simplify algebra.

### Several missing indices

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

### Empty broadcasts

| Input | Required result |
| ----- | --------------- |
| `index_sum(indexed(broadcast_to(square_root(Integer(-1)), (0,)), (i,)), (i,))` | Value zero with shape `()` and no free indices. Neither the empty broadcast nor the empty sum evaluates the real square root. The prototype remains subject to static construction checks. |

This does not make `square_root(Integer(-1))` a defined real value. Symbolic
traversal under SH19 is distinct from evaluating an array member.

## Maintaining rule references

Rule IDs identify rewrite and shared-requirement rows, including retained cases.
They stay with their rules across moves and are never renumbered or reused.
New rules receive new IDs; prefixes identify families. Contracts and schedules
are unnumbered. Shared conditions and notation apply to every row.
See the [renumbering record](../scratch/language/rule-renumbering.md) for former IDs.

ST05 (retained nested positional application) is retired: full `Apply` has no
remaining argument slots.
