# Rewriting and lowering

Rewriting replaces expressions without changing their mathematical value.
Lowering applies rewrites to meet the requirements of later computations.

## Choosing a representation

A *contract* describes a required representation, such as expanded derivatives.
Rules give individual rewrites; recipes compose them to meet contracts.

`lower` is one recipe for reaching the [core normal form](#the-core-normal-form).
Other recipes may reach that form or satisfy different contracts, subject to
the [shared requirements](rewriting.md#shared-requirements).
Core normal form leaves reference-cell integrals, terminal gradients and function
evaluations for a library or compiler to handle.

| Page | Contents |
| --- | --- |
| [Rewrite machinery](rewriting.md) | Common traversal, rule selection, shared linearity rules and type completion |
| [Substitution](substitution.md) | Binding scopes, capture avoidance and fresh names |
| [Differentiation](differentiation.md) | Directional and spatial derivative recipes, retained gradients and derivative formulas |
| [Applications and mappings](lowering.md) | Applications, mesh integration, pull-backs and tensor arithmetic |
| [Expression optimisation](optimisation.md) | Optional simplification, numerical contract and examples by node class |

## Standard lowering recipe

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `lower(f)` | | `f: AbstractExpression` | An equivalent expression in core normal form |

| `lower` condition | Requirement or result |
| ----------------- | --------------------- |
| Input | Core expressions and [prepared geometry](rewriting.md#evaluation-traversal); extension nodes need an extension recipe first |
| Unsupported rewrite | `TypeError` |
| Stable output violates core normal form | `ValueError`, identifying offending nodes |
| A complete spatial round returns to an earlier state | `ValueError`; a cycle is not completion |

Run these passes in order. Derivative preparation may expand integrals and
spatial mappings locally (D21–D23). Spatial passes repeat because pull-backs
can expose new gradients.

| Recipe step | Required result |
| ----------- | --------------- |
| [`expand_applications`](lowering.md#expanding-applications) | Applications expanded |
| [`apply_derivatives`](differentiation.md#applying-derivatives), including application cleanup | Applications and derivatives expanded |
| [`expand_integrals`](lowering.md#expanding-integrals) | Also mesh integrals expanded |
| Repeat [`apply_gradients`](differentiation.md#applying-gradients), then [`expand_pull_backs`](lowering.md#expanding-pull-backs), with application cleanup and reconsideration under SH19–SH24 | Stop under the conditions below |
| [`expand_arithmetic`](lowering.md#expanding-tensor-arithmetic) | Core normal form |

Compare states after complete spatial rounds using structural equality,
including bound names and labels. Fresh naming alone is not progress.

| Round result | Next step |
| ------------ | --------- |
| Neither spatial pass nor application cleanup changed the expression; all contracts except tensor arithmetic expanded hold | Proceed to `expand_arithmetic` |
| Neither changed it; another core normal-form contract fails | Raise `ValueError` with the remaining offending nodes |
| Some step changed it, but the final state equals an earlier round boundary | Report a cycle; do not treat it as an unchanged round |
| Otherwise | Run the next round |

## Representation checks

### Representation contracts

| Contract | Required representation |
| -------- | ----------------------- |
| Core expressions | Structurally valid concrete core nodes and downstream terminals |
| Applications expanded | Every `Apply` applies a `Function`, downstream terminal, `Gradient` or `PullBack` |
| Derivatives expanded | No `DirectionalDerivative` |
| Mesh integrals expanded | No `Integral` over a `Mesh` |
| Gradients expanded | Every `Gradient` matches [Retained gradients](differentiation.md#retained-gradients); every `Divergence` or `Curl` matches [Retained Piola patterns](differentiation.md#retained-piola-patterns) under no cell pull-back |
| Pull-backs expanded | Every `PullBack` matches [Retained pull-backs](lowering.md#retained-pull-backs) |
| Tensor arithmetic expanded | Every `Sum`, `Product` and `Division` has scalar value shape |
| Reference integration | Derivatives expanded and mesh integrals expanded |
| Core normal form | The allowed patterns in [The core normal form](#the-core-normal-form) |

Checks report offending nodes and reasons without rewriting. Preservation
requirements cover domains and operators; quadrature and tabulation have
separate contracts.

### The core normal form

Core normal form combines all the expansion contracts above with core
expressions. Its allowed nodes are:

| Chapter | Node class | Allowed in the core normal form |
| ------- | ---------- | ------------------------------- |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Integer`, `Float`, `Complex`, `Zero` | Always |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Sum`, `Product`, `Division` | Scalar value shape |
| [Numbers and scalar operations](expressions.md#numbers-and-scalar-operations) | `Power`, `Absolute`, `Conjugate`, `RealPart`, `ImaginaryPart`, `SquareRoot`, `Logarithm` | Always |
| [Tensors and indexing](expressions.md#tensors-and-indexing) | `Tensor`, `BroadcastTo`, `Indexed`, `ComponentTensor`, `IndexSum` | Always; index simplification is left to the compiler |
| [Tensor operations](expressions.md#tensor-operations) | `Inverse`, `Determinant` | Always, for the compiler to evaluate |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Coordinate` | Always |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Gradient` | A [retained gradient](differentiation.md#retained-gradients) |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Divergence`, `Curl` | A [retained Piola pattern](differentiation.md#retained-piola-patterns), under no pull-back by a cell |
| [Calculus on a domain](geometry.md#calculus-on-a-domain) | `Integral` | Over a domain that is not a `Mesh`, such as a reference cell or `EuclideanSpace(0)` |
| [Maps between domains](topology.md#maps-between-domains) | `PullBack` | A [retained pull-back](lowering.md#retained-pull-backs) |
| [Fields on a mesh](geometry.md#fields-on-a-mesh) | `MappedNormal` | Always |
| [Fields on a mesh](geometry.md#fields-on-a-mesh) | `PushForward` | Under no pull-back by a cell |
| [Function spaces](forms.md#function-spaces) | `Function` | Always |
| [Variational forms](forms.md#variational-forms) | `Form` | Always, including a form whose expression is a form (see [`Form`](forms.md#form)) |
| [Variational forms](forms.md#variational-forms) | `Apply` | Its form is a `Function`, a downstream terminal, a `Gradient`, or a `PullBack`, with each operand in the core normal form |
| | A terminal of a downstream library | Always |

The normal-form check follows the [evaluation traversal](rewriting.md#evaluation-traversal);
it does not check stored geometry or labels.

Core normal form does not imply algebraic simplification or an efficient
evaluation order. It may retain component-tensor selections, identity-tensor
contractions, redundant computations and unfactorised sums. Consumers may
simplify or reorganise these expressions under the
[optimisation contract](optimisation.md); no particular
optimisation or schedule is required.

[Compiling numerical kernels](compiler.md) describes possible consumer choices.

Assembly code receives the form layer (see [Grammars](rewriting.md#grammars)), including
sums and products of forms, real and complex zeros, nested forms, and applications
on supplied or spatially transformed forms. The static type gives the argument
order (see [`Form`](forms.md#form)).

A form compiler receives the integrand layer, consisting of integrals over reference cells
and points, with arithmetic and index notation over constants, reference
coordinates and the library reads of [Functions in a space](forms.md#functions-in-a-space).

```{dropdown} Reusing geometry calculations
Jacobians, inverse Jacobians, densities and side normals remain subexpressions.
A compiler can compute them once per set of map data and reuse them across integrals.
```

## Progress and stability

### Progress of the reference recipe

Expression size may grow, so node count is not a progress measure. The rules
remove an operation, move it to operands, or stop at an explicit retained case.

| Rule family | Progress |
| ----------- | -------- |
| Applications and directional derivatives | Apply or differentiate the current node, then process operations on its operands. An inner derivative is processed before an outer one. Retained applications are not expanded again. |
| Mesh integration | Replace a mesh integral by reference-cell integrals. Repeating the pass cannot expand those as mesh integrals. |
| Gradient reconstruction | Apply ST08 only to a gradient outside ST12. Component rules descend through operands. ST09 introduces a gradient that matches ST12, so reconstruction cannot wrap it again. |
| Application transport | ST07, ST10 and ST11 move applications outside spatial operations. ST03 prevents application expansion from reversing this move. |
| Pull-back algebra | Move a pull-back into operands, compose an existing map, or retain a read. Composition exposes the existing map's parameterisation for processing; it does not recreate the original map composition. |
| Cell transport and Piola identities | Replace the matched cell pull-back by its transformation formula. New gradients are handled in the next spatial round and stop at ST12 where appropriate. |
| Tensor arithmetic | ST13–ST16 replace a tensor arithmetic node by scalar arithmetic and index notation. They introduce no tensor arithmetic nodes. |
| Retained cases | Do not create fresh indices or schedule another rewrite of the unchanged node |

These steps use the finite, acyclic input and stored geometry. Other recipes
must justify their own progress; cycle checks alone do not prove termination.

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

[Geometry differentiation](differentiation.md#geometry-preparation) may require spatial
preparation before differentiation and cleanup afterwards. Index simplification,
such as reducing an `Indexed` of a `ComponentTensor`, belongs to the compiler.
