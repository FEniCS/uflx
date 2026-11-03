# Extensions

## Extending the core

The core has a small, fixed set of node classes and an open set of passes. A
downstream library extends UFLx mainly by adding passes, not by adding
properties to the node classes of the core. A new pass needs one rule for each
node class of the core and does not change the core. Only a language extension
adds node classes.

The core supports five kinds of extension:

- [New questions](#new-questions), such as physical units, polynomial degree
  or linearity checks.
- [New subclasses](#new-subclasses), such as the elements and functions of
  Basix, or new coordinate domains.
- [New operations](#new-operations), such as the trace or the cross product.
- [New transformations](#new-transformations), such as compiler lowering of
  index notation.
- [Language extensions](#language-extensions), such as external operators.

A property that the rules of the core depend on is part of the core, not an
extension. For example, the scalar type is part of the core, because
`Absolute` is real for a complex operand. Adding such a property changes the
core.

```{dropdown} Legacy UFL
In legacy UFL, properties are stored on the node classes, and the set of properties is fixed. Every new property touches every node class, so features such as complex numbers and dual spaces were patched in across the node classes.
```

### Language extensions

A language extension is a package outside the core that adds non-leaf node
classes. Each extends `AbstractExpression`, never a concrete core class
(see [Leaves and subclassing](rules.md#leaves-and-subclassing)), and provides:

- the value of the node, as the section of a lazy node class of the core
  states the value (see [Evaluation](evaluation.md#evaluation)), and
- its expression operands for the [evaluation traversal](evaluation.md#evaluation-traversal),
  distinguishing them from metadata, and
- a lowering to a mathematically equal expression of core node classes, if
  such an expression exists.

An extension node may state argument spaces, making it an operator with
arguments. The extension's action-evaluation table defines the action on
each argument; dual arguments alone do not require an extension node.

#### Rule tables

A pass is built from an ordered `tuple` of rule tables, each a list of rules
pairing a rule pattern with a rewrite. A rule pattern may match the node
class, operand classes, a condition on the static type, or the pass result
on an operand. The first matching rule wins, in table order and then rule
order. Each rule receives the whole pass, including every table, to apply to
operands.

- A table of a language extension may match only some nodes. The table of
  the core comes last, and has a rule for every node class of the core.
- A pass applied to a node that no rule matches raises `TypeError`.
- The lowering table of a node class with a lowering matches the nodes of the
  class, lowers the node and applies the pass to the lowering.
- The caller chooses the tables and their order, including extension tables;
  there is no global registry. The result depends on this order, not on the
  order in which extensions are loaded.

An extension supplies its own rules or its lowering table for each supported
pass, whether from the core or another library, such as `units`. Earlier
extension tables may override core rules. Rewrite rules preserve value.
Every rule must agree, up to mathematical equality, with lowering followed
by the pass.
An identity of a composition, such as the identity of Piola, needs node
classes for that composition, as with `Divergence` and `Curl` in the core
(see [`PushForward`](geometry.md#pushforward)).

#### Examples

- A language extension for local external operators adds a node class for an
  operation known only by name, such as `exp`, a material law or a neural
  network. The node class has no lowering. The table of the extension for a
  derivative gives the derivative of the operation, and the table for a
  pull-back gives the node of the pulled-back operands, because the operation
  is pointwise. A pass with no table of the extension raises `TypeError` for
  the node.
- An interpolation node may retain the target space or delay expansion into
  sampling and basis reconstruction. Its lowering is the library operation
  of [Cell-local interpolation](#cell-local-interpolation). Forward
  interpolation itself needs no new node class.
- A named composition that must survive a pass, or print as written, is a node
  class of a language extension whose lowering is the composition. A
  `Transpose` node class lowers to `component_tensor(A[k, l], (l, k))`. The
  table of the extension for a pull-back gives the `Transpose` of the
  pull-back, because the transpose is pointwise, so `Transpose` survives the
  pull-back. A pass such as a derivative takes the lowering table of
  `Transpose`.
- For a pass that reduces an `Indexed` of a `ComponentTensor`, the table of
  the `Transpose` extension matches an `Indexed` of a `Transpose`, and
  rewrites `Transpose(A)[i, j]` as `A[j, i]`. The table comes before the table
  of the core, which has a rule for every `Indexed`.

[Shape derivatives](#shape-derivatives) give an example of an extension whose
pass runs between two lowering stages.

```{dropdown} Open: operator composition for preconditioning
A linear preconditioner maps residuals in $V^\star$ to corrections in $V$.
An auxiliary-space construction is $P B_W P^\star$, with $P: W \to V$ and
$B_W: W^\star \to W$. Dual arguments and typed operator composition could
express this construction. It is open whether arbitrary operator-valued
arguments are also useful, or composition belongs in the solver layer.

A concrete finite element implementation should establish how the symbolic
operations connect to assembly, matrix-free application and auxiliary solves
before extending `Form` further. Preconditioner construction need not be
linear in its operator inputs, even when its residual-to-correction action is
linear.
```

### New questions

A downstream library asks a new question of an expression with a pass. The pass
has a rule table for the node classes of the core and the terminal classes of the
library. A pass that meets a node class of a language extension also has a
table for that class before the table of the core, such as the lowering table
(see [Rule tables](#rule-tables)).

For example, a units library assigns metres to `X = coordinate(RD(2))`
and treats unlisted terminals as dimensionless. Its pass derives square
metres for `X[0] * X[1]`, a dimensionless value for `gradient(X[0], RD(2))`,
and an error for `X[0] + X[0] * X[1]`. The library calls the pass when needed.

### New operations

A downstream library builds a new operation from existing node classes, as a
construction free function. As in the core, it may read the static type of
the operands, such as a value shape, to choose the composition (see
[Node class or free function](expressions.md#node-class-or-free-function)). For example,
`cross(a, b)` checks that both operands have value shape `(3,)` and returns:

```python
tensor(
    (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )
)
```

`trace(A)` reads no extent,
because `index_sum(A[k, k], (k,))`, with a derived index `k`, is the trace for
every extent. Operations such as `outer(u, v)` read operand ranks to build
index tuples, without inspecting axis extents (see
[Tensor algebra](../examples/tensors.md#tensor-algebra)).

#### Cell-local interpolation

`interpolate(f, V)` is a library construction free function. Its value is the
field in `V` whose degrees of freedom are given by applying the target's
sampling rule to `f`.

| `interpolate(f, V)` | |
| ------------------- | --- |
| Takes | `f: AbstractExpression`, `V: AbstractFunctionSpace` |
| Requires | `V` a finite element space on a reference cell, `f` with that domain or no domain and the value shape and scalar type of `V`, both with argument spaces `()`, and `f` defined at every sampling point |
| Static type | The domain, value shape and scalar type of `V`, the free indices of `f`, and argument spaces `()` |

The library of `V` defines each degree of freedom as a weighted combination
of [point evaluations](geometry.md#point-evaluation), including any quadrature for
moments. The supplied rule defines interpolation even when the moments are
not exact. Interpolation is linear for a fixed rule and geometry.

For a target obtained from a reference element by a
[transformation law](geometry.md#transformationlaw), interpolation on a physical cell
applies the target's inverse law to the pulled-back input, then maps the
reference interpolant to the physical cell by the target's law.
Geometry-dependent basis construction and interpolation belong to the library
and are outside this description.
A `PullBack` alone does not transform value components.

See [Dual transfer through interpolation](../examples/forms.md#dual-transfer-through-interpolation)
for a pairing with a target load.

Only value sampling and affine, invertible cell maps between equal dimensions
are covered. Derivative sampling is not yet covered.

```{dropdown} Open: an interpolation node
Whether a language extension node should retain membership of `V` or delay
expansion is open.
```

### New subclasses

An abstract class states the contract that downstream subclasses fulfil
without changing the core. Passes handle them through that contract, which
belongs to the class, never a global registry. For example, a downstream
`AbstractCoordinateDomain`, such as
[`UnitSimplex`](../examples/library.md#the-unit-simplex), provides the dimension; passes
handle its coordinate function through `Coordinate`.

Information that a downstream library attaches to the graph goes in one of two
places:

- Information that defines the variational problem goes in the construction
  arguments of a leaf. For example, the space that a function belongs to
  changes the meaning of an expression, so the space is part of the leaf.
- Information that only checks or analyses the variational problem may stay in
  a table of the library, keyed by leaves. For example, units check an
  expression but do not change the meaning of the expression.

### New transformations

A new transformation is a pass with one rule for each node class of the core
and one rule for each terminal class that the transformation rewrites. For
example, a compiler may reduce an `Indexed` of a `ComponentTensor`, following
the [substitution and binding](evaluation.md#substitution-and-binding) rules.

### Shape derivatives

> A language extension, not a rule of the core.

`ShapeDerivative(f, x, h)` is the derivative of $f$ when the coordinate function
$x$ of a mesh `image(x, Q)` moves in the direction $h$ (see [A shape
derivative](../examples/surfaces.md#a-shape-derivative)). The pass of the extension runs
between two spatial-lowering stages, where $x$ occurs explicitly in the
parameterisations of the cells, and lowers the node to `DirectionalDerivative(f,
x, x, h)`, which holds the domains fixed. The extension recipe first exposes
geometry while retaining its own nodes, applies shape differentiation, then
re-enters `lower` to finish derivatives, resulting actions and spatial operations
and check core normal form. The first stage is not a completed call to core
`lower`, since extension nodes remain. A function on the topology keeps its
values on the reference cell, so the function moves with the mesh.

The direction is any direction that the core derivatives allow, an expression
on the topology of value shape `(n,)` or a constant. With an argument as the
direction, the Gateaux derivative gives the shape gradient. An ALE mesh is
`image(x0 + d, Q)`.
