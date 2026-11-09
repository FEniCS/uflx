# Geometry and calculus

A local cell configuration acquires physical geometry through a coordinate expression:

```python
# T describes one triangle slot; the package supplies its local data.
Q = RD(2)                       # Ambient space: the plane.
X = Lagrange(T, 1, (2,))         # Library-supplied space of coordinate functions.
x = function(X, "x")            # Physical positions, with values supplied by the solver.
D = mesh(x, Q)                  # Physical mesh; its topology comes from x.
```

Changing the coordinate expression can change the physical geometry while
retaining the topology and the functions defined on it. Fields on the topology
are carried to the physical mesh by [push-forward](#pushforward).

The only metric is the ambient Euclidean metric, which defines the tangential
[gradient](#gradient), integration [density](#jacobian-and-density) and unit
[normal](#normal).

## Meshes

### `Mesh`

`Mesh(x, Q)` is one physical local field domain, with coordinate expression $x$ on
an `AbstractTopology` $T$ and ambient domain $Q$. Cells keep their identities
when their physical coordinates coincide. The mesh stores no
global connectivity or coordinate values. The package supplies the local data.
Its finitely many cell slots describe one kernel configuration, not a runtime mesh.

| `Mesh(x, Q)` | |
| ------------ | - |
| Extends | `AbstractDomain` |
| Free function | `mesh(x, Q)` |
| Takes | `x: AbstractExpression`, `Q: AbstractCoordinateDomain` |
| Properties | `coordinate_function` $x$, and `ambient_domain` $Q$ |
| Requires | $x$ with domain $T$, a topology whose cell types are reference cells; $x$ real, with value shape `(n,)` for the ambient coordinate dimension of $Q$, free indices `()` and argument spaces `()`; `T.dimension <= n` |
| Derived | `topology` $T$ = `x.domain`, `dimension` = `T.dimension`, and `pieces`, below |
| Raises | `TypeError` if the coordinate expression has no topology domain, and `ValueError` if another requirement fails |

`mesh(x, Q)` returns `Mesh(x, Q)`. The topology is determined by $x$ and is
not a separate construction argument. Mesh identity depends on $x$ and $Q$;
replacing the coordinates reconstructs the mesh from the new coordinates.

A mesh describes one geometry, not an indexed family of geometries. Coordinates
may depend on symbolic coefficients, deformation fields or scalar parameters,
but have no free indices. Construction and reconstruction enforce this rule.
A solver may evaluate a form for different coordinate values or manage several
meshes outside the language. Shape differentiation does not require an indexed
mesh (see [Geometry dependence](#geometry-dependence)).

Deformed surfaces, isogeometric patches and analytic maps are meshes of
coordinate functions. Continuity across cells is a property of the space of
$x$. Independent patches or points that support fields use a topology with
independent cells; no sides are required. An indexed family used only for
integration or sampling needs a `DomainMap`, not a mesh (see
[Mapped integration families](#mapped-integration-families)).

#### Cell maps and pieces

For slot $j$ of $T$, let $K_j$ be its reference cell. `pieces` is the derived
`tuple` of maps
`domain_map(K_j, pull_back(x, cell_map(T, j)), Q)`.
These maps are views of the local mesh, not construction arguments. They have
no free indices and introduce no extra geometry parameters.

The ambient domain $Q$ is a coordinate domain. A mesh embedded in another
mesh has its own physical coordinates; cell or side maps specify its
relation to the other mesh's cells.

#### Mesh cell maps

`cell_map(D, j)` returns `CellMap(D, j)`, whose target is the mesh $D$.
Its parameterisation is `pull_back(x, cell_map(T, j))`. This expression
is derived, not a construction argument. Substitution reconstructs the map
from $D$ and $j$, preserving the selected cell independently of its
coordinates. Piece $j$ instead
has target $Q$ and records only physical coordinates. These maps have the
same parameterisation but different meanings: only the cell map selects a
cell of $D$ for evaluating a mesh field.

For a field $f$ on $D$, `pull_back(f, cell_map(D, j))` is its local
expression on $K_j$ for local slot $j$. Pulling back
`coordinate(D)` gives `pull_back(x, cell_map(T, j))`.

#### Examples

Use the two-triangle topology of the [static topology
example](../examples/library.md#a-static-topology), on reference cells:

```python
K = ReferenceSimplex(2)
square = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)), ())
```

Let `x` be a vector Lagrange function of degree 1 on `square`, of value shape
`(2,)`, with vertex values $(0, 0)$, $(1, 0)$, $(0, 1)$ and $(1, 1)$.
The resulting physical mesh is the square:

```python
D = mesh(x, RD(2))
```

- `D.topology` is `square`, with two triangle slots. `cell_map(D, 0)`
  reads the first; `cell_map(D, 1)` reads the second.
- `pull_back(coordinate(D), cell_map(D, 1))` is
  `pull_back(x, cell_map(square, 1))`, which the library evaluates to
  $(1 - X_1, X_0 + X_1)$ for the supplied coordinates.

For a local configuration `M` with types `(Kt, Kq)`, `mesh(x, RD(2))` has
one triangle and one quadrilateral piece. A field `u` on `M` has local reads:

```python
u_tri = pull_back(u, cell_map(M, 0))
u_quad = pull_back(u, cell_map(M, 1))
```

The package can reuse this configuration for triangle–quadrilateral
interfaces. Separate cell kernels use single-slot configurations. A library
specialises forms for the configurations needed by a mixed mesh.

```{dropdown} Legacy UFL
In legacy UFL, a `Mesh` is `mesh(x, Q)` for the coordinate element of `x`, and `SpatialCoordinate(mesh)` is `coordinate(D)`.
```

## Calculus on a domain

### `Coordinate`

`Coordinate` maps each point of a domain to its coordinates in the ambient
coordinate domain. It is a terminal.

| `Coordinate(P)` |                                                    |
| --------------- | -------------------------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `coordinate(P)`                                    |
| Takes           | `P: AbstractDomain` |
| Properties      | `domain`, an `AbstractDomain` $P$                  |
| Static type     | Value shape `(n,)`, real, the domain $P$          |

$n$ is the dimension of the ambient coordinate domain of $P$. For a
coordinate domain $P$, the coordinate function is the identity,
$X : P \to R^n$, $X(p) = p$, and the gradient of $X$ is the identity matrix
(see [`Gradient`](#gradient)).

For a mesh $D$, `coordinate(D)` is the physical coordinate field on $D$.
Its pull-back by a cell map is that cell's parameterisation (see
[`Mesh`](#mesh)). This differs from `D.coordinate_function`, which is
the stored coordinate expression on the topology.
For a mesh `mesh(x, Q)`, the pull-back by `cell_map(D, j)` is
`pull_back(x, cell_map(T, j))`.

For example, `coordinate(RD(3))` has value shape `(3,)` and domain `RD(3)`.

### Valid pair

An expression and a domain $(f, P)$ are a valid pair if $P$ is an
[`AbstractDomain`](topology.md#abstractdomain) and $f$ has the domain $P$ or no domain.
The type of $P$ is checked statically; construction raises `ValueError` if
the domains do not match. For $f$ with no domain, use the
[constant-field interpretation](expressions.md#domain-free-expressions) on $P$.

### `Gradient`

`Gradient` is the gradient of an expression on a domain $P$. The gradient adds a
last axis in the coordinate frame of the ambient coordinate domain of $P$, of
dimension $n$.

| `Gradient(f, P)` |                                                                                                   |
| ---------------- | ------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `gradient(f, P)`                                                                                  |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `domain` $P$, with $n$ the dimension of the ambient coordinate domain of $P$ |
| Requires         | A [valid pair](#valid-pair) $(f, P)$                                                              |
| Static type      | Value shape `f.value_shape + (n,)`                                                                |
| Raises           | `ValueError` for $f$ if a requirement fails                               |

The gradient acts on each component of an array-valued $f$. Entry `j` of the
last axis of `Gradient(f, P)` is the partial derivative $\partial_j f$ in the
coordinate $X_j$ of the frame.

On a mesh, `Gradient`, `Divergence` and `Curl` act cellwise.
They do not include distributional interface terms for discontinuous fields;
traces and interface integrals are written explicitly.

#### The frame

For a domain $P$ of dimension $k$, the frame is the coordinate frame of the
ambient coordinate domain of $P$, of dimension $n$. A coordinate domain has
$k = n$.

- For $k = n$, $\partial_j$ is the partial derivative in the coordinate
  $X_j$ of $P$.
- For $k < n$, such as a mesh that is a surface in $R^3$, the
  coordinate directions of the ambient coordinate domain are not tangent to
  $P$. `Gradient(f, P)` is then the tangential gradient, the gradient of an
  extension of $f$ off $P$, projected orthogonally onto the tangent space of
  $P$ with the Euclidean metric of the ambient coordinate domain. The
  tangential gradient does not depend on the extension.

The divergence and the curl are [`Divergence` and `Curl`](#divergence-and-curl),
defined by `Gradient`. Other operators of vector calculus are fixed
compositions of `Gradient`, which a library provides as free functions (see
[New operations](extensions.md#new-operations)).

#### Pull-backs of `Gradient`

The gradient depends on the coordinates of $P$, so a pull-back does not
commute with `Gradient`. `PullBack(Gradient(f, Q), D)` is the gradient of $f$
on $Q$, evaluated at $\psi(p)$.

#### Examples

Start with a scalar expression on `RD(3)`:

```python
P = RD(3)
X = coordinate(P)
x, y, z = X[0], X[1], X[2]
f = x * y * z
```

| Expression                              | Value shape | Meaning                                       |
| --------------------------------------- | ----------- | --------------------------------------------- |
| `gradient(f, P)`                        | `(3,)`      | The gradient $(yz, xz, xy)$, after evaluation |
| `gradient(gradient(f, P), P)`           | `(3, 3)`    | The Hessian of $f$                            |
| `gradient(X, P)`                        | `(3, 3)`    | The identity matrix, after evaluation         |

- Let `M` be a mesh of dimension 2 in `RD(3)`, and `w` a field
  on `M`. `gradient(w, M)` is the surface gradient of `w`, with value shape
  `(3,)`. `gradient(coordinate(M), M)` is the orthogonal projection onto the
  tangent plane of `M`, with value shape `(3, 3)`.
- `gradient(X, RD(2))` raises `ValueError`, because `X` has the domain
  `RD(3)`.

### `Divergence` and `Curl`

`Divergence` and `Curl` are the divergence and curl of an expression $u$ on a
domain $P$, with $n$ the dimension of its ambient coordinate domain.

| `Divergence(u, P)` |                                                                     |
| ------------------ | ------------------------------------------------------------------- |
| Extends            | `AbstractExpression` |
| Free function      | `divergence(u, P)`                                                  |
| Takes              | `u: AbstractExpression`, `P: AbstractDomain` |
| Properties         | `expression` $u$ and `domain` $P$                                   |
| Requires           | A [valid pair](#valid-pair) $(u, P)$, and $u$ of rank at least 1 with the last extent $n$ |
| Static type        | Value shape `u.value_shape[:-1]`                                    |
| Raises             | `ValueError` for $u$ if a requirement fails |

| `Curl(u, P)`  |                                                                     |
| ------------- | ------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `curl(u, P)`                                                        |
| Takes         | `u: AbstractExpression`, `P: AbstractDomain` |
| Properties    | `expression` $u$ and `domain` $P$                                   |
| Requires      | A [valid pair](#valid-pair) $(u, P)$, $n$ 2 or 3, and $u$ of value shape `(n,)` |
| Static type   | Value shape `(3,)` for $n = 3$ and `()` for $n = 2$                 |
| Raises        | `ValueError` for $u$ if a requirement fails |

The value of `Divergence(u, P)` is the trace of `gradient(u, P)` over the last
axis of $u$ and the axis of the gradient. With `G = gradient(u, P)`, so
$G_{ij}=\partial_j u_i$, `Curl(u, P)` is the antisymmetric contraction:

| Ambient dimension | Value |
| ----------------- | ----- |
| 2 | `G[1, 0] - G[0, 1]` |
| 3 | `tensor((G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]))` |

`Gradient` stays the one derivative that defines the others. The two node
classes are in the core so that the identity of Piola
holds structurally (see [`PushForward`](#pushforward)).

### `DirectionalDerivative`

`DirectionalDerivative` differentiates an expression $f$ with respect to a
variable $x$ at a point $a$ in a direction $h$,

$$
\mathrm{D}(f, x, a, h) = \frac{d}{d\varepsilon} f[x := a + \varepsilon h] \Big|_{\varepsilon = 0},
$$

with $\varepsilon$ real. Substitute free occurrences of $x$ in a compatible
fixed-space representation of $f$ with [geometry dependence](#geometry-dependence) exposed.
For tuple variables, matching tuples $a$ and $h$ supply all replacements
$x_k:=a_k+\varepsilon h_k$ simultaneously. Members may mix real and complex
spaces. The functions of $x$ are bound in $f$, not in $a$ or $h$.

| `DirectionalDerivative(f, x, a, h)` |                                                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Extends                          | `AbstractExpression` |
| Takes                            | `f: AbstractExpression`, `x: Function \| tuple[Function, ...]`, `a, h: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties                       | `expression` $f$, `variable` $x$, `point` $a$ and `direction` $h$ |
| Requires                         | $x$ a [bindable function](forms.md#bindable-functions), or a `tuple` of at least two distinct bindable functions; member by member, $a$ and $h$ admissible values for the space of $x$ as in [argument admissibility](forms.md#admissible-argument-values); variable spaces independent of all variable members |
| Static type                      | The value shape of $f$, the free indices of $f$ followed by those of the members of $a$ and $h$, each `Index` taken once, and the argument spaces and the scalar type of $f$ |
| Raises                           | `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `directional_derivative(f, x, a, h)` | | `f: AbstractExpression`, `x: Function \| tuple`, `a, h: AbstractExpression \| tuple` | `DirectionalDerivative(f, x, a, h)` |

A real-linear differential is assumed, even for complex variables and
directions. Real points and directions enter complex spaces by
[complex embedding](evaluation.md#type-compatibility-of-rewrites).
Dual variables accept symbolic forms with matching functional signatures.

This is the sole derivative primitive for function variables. Binding and
complex projection give [derivative forms](forms.md#derivative-forms).
A domainless direction has zero spatial gradient, so varying `u` along it
leaves `gradient(u, P)` unchanged when `P` is independent of `u`.

For `u`, `a` and `h` in a real scalar space of constants,
`directional_derivative(u * u, u, a, h)` has value $2 a h$.

An indexed direction gives a family of derivatives. For scalar constants
`h0` and `h1` in the same space:

```python
i = index("i")
directions = tensor((h0, h1))
d = directional_derivative(u * u, u, a, directions[i])
```

`d` has free index `i` and value `Integer(2) * a * directions[i]`. Tensor derivatives
use families of unit directions; the directions specify any symmetry. See
[Hyperelasticity and the stress](../examples/forms.md#hyperelasticity-and-the-stress).

#### Geometry dependence

Differentiation counts dependencies through domains, maps, spaces and library
reconstruction, holding independent field data fixed. Before applying ordinary
chain rules, expose dependent geometry on fixed reference domains or use an
explicit library derivative rule. Independent geometry needs no lowering.
See [Geometry preparation](evaluation.md#geometry-preparation).

Variable spaces must be independent of the variables. The differentiated
expression requires a fixed input signature and output domain compatible with
its declared type. A valid form whose argument spaces depend on free geometry
must first be expressed on fixed data spaces. A field on a moving domain
requires an explicit identification, such as a pull-back to a reference domain.
Differentiation chooses no implicit transport or interpolation.

Libraries must expose geometry-dependent reconstruction as symbolic operands
or derivative rules, including basis and coefficient transformations. Hidden
runtime dependencies violate this contract. Missing reconstruction or derivative
rules raise `ValueError` identifying the unresolved dependency, never a silent
zero derivative.

Topology, connectivity, index extents and discrete element choices remain fixed.
Derivatives assume admissible, differentiable geometry; topology changes and
singular or nonsmooth configurations are excluded.

#### Nested differentiation

Resolve inner derivatives first, including their geometry dependence. Their
points and directions lie outside the variable scope, so outer derivatives
include those operands' free dependencies. Fixed-direction Hessians use
independent direction functions. Mixed derivatives commute only under the
required smoothness and fixed-direction hypotheses.

### `Integral`

`Integral` integrates an expression $f$ over a domain $P$. A basic integral
needs no topological notion of the domain.

| `Integral(f, P)` |                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `integral(f, P)`                                                                            |
| Takes            | `f: AbstractExpression`, `P: AbstractDomain` |
| Properties       | `expression` $f$ and `integration_domain` $P$                                                           |
| Requires         | A [valid pair](#valid-pair) $(f, P)$, or $P$ a mesh and $f$ with the domain `P.ambient_domain` |
| Static type      | The value shape and free indices of $f$, with no domain |
| Raises           | `ValueError` for $f$ if a requirement fails                         |

For a coordinate domain $P$ of dimension $n$, `Integral(f, P)` is the
Lebesgue integral in the coordinates $X$ of $P$,

$$
\int_P f \,\mathrm{d}X.
$$

For $n = 0$, the integral is the value of $f$ at the single point of $R^0$.

On a mesh, the integral sums the integrals over its pieces with
their induced densities (see [the expansion](#integral-over-a-mesh)).

An array-valued $f$ integrates component by component. A later pass might
approximate an `Integral` over a coordinate domain by
[quadrature](#quadrature), or evaluate it exactly when the rule is exact.

#### Examples

With `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0], T)` is $\int_T X_0 \,\mathrm{d}X$. A later pass evaluates
  the integral to $1/6$.
- `integral(Integer(1), T)` is the area $1/2$ of `T`.
- `integral(X[0], T) * integral(X[1], T)` is a product of two integrals, with
  value shape `()` and no domain.

## Maps between domains

### `DomainMap`

`DomainMap` is a map $\psi$ from a parameter domain $P$ into an ambient coordinate domain
$Q$, with geometric image $\psi(P)$. The sections below write $D$ for the map
or its parameterised image, retaining multiplicity. A `DomainMap` is not an
`AbstractDomain` and is not the domain of an expression.

| `DomainMap(P, psi, Q)` |                                                                                                                                                                                                             |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                          | `Base` |
| Free function                    | `domain_map(P, psi, Q)`                                                                                                                                                                           |
| Takes                            | `P: AbstractCoordinateDomain \| AbstractTopology`, `psi: AbstractExpression`, `Q: AbstractCoordinateDomain` |
| Properties                       | `parameter_domain` $P$, an [`AbstractCoordinateDomain`](topology.md#abstractcoordinatedomain) or an [`AbstractTopology`](topology.md#abstracttopology) of dimension $k$, `parameterisation` $\psi$, and `ambient_domain` $Q$, an `AbstractCoordinateDomain`                                                   |
| Requires                         | $\psi$ with the domain $P$ or no domain, $\psi$ real with value shape `(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$ and the argument spaces `()`, and $k \le n$ |
| Derived                          | `free_indices`, the free indices of $\psi$ |
| Raises                           | `ValueError` if another requirement fails                                                                                                                                   |

The topological dimension is `P.dimension`; the geometric dimension is
that of the ambient coordinate domain of $Q$. Maps need not be injective.
Coordinates alone do not select a cell or reference preimage in a mesh;
mesh-field evaluation uses [cell or side maps](#a-family-of-points-in-a-mesh).

A $\psi$ with free indices makes $D$ a family of domains, with one domain for
each value of the free indices. For example, with an array `x` of value shape
`(m, n)` and an index `i`, `domain_map(RD(0), x[i], RD(n))` is a family of $m$
points in `RD(n)`, retaining distinct members even at identical coordinates.

For example, the hypotenuse of the unit triangle has parameter domain
$P = [0,1]$ and parameterisation $\psi(t) = (1-t,t)$. Choose `Q = RD(2)`
for the edge in the plane, or `Q = UnitSimplex(2)` for the edge in the triangle.

### `PullBack`

`PullBack` is the pull-back $\psi^* f$ of an expression $f$ on $Q$ by a map $D$
with parameter domain $P$ and ambient domain $Q$. A `DomainMap` additionally
specifies its coordinate parameterisation $\psi$.

Here $P$ and $Q$ denote the domains from which and into which the map goes.
Their API properties depend on the map's role:

| Map class | Starting domain $P$ | Target $Q$ |
| --------- | ------------------- | ---------- |
| `DomainMap` | `parameter_domain` | `ambient_domain` |
| `CellMap` | `reference_cell` | `parent` |
| `SideMap` | `restriction_domain` | `parent` |

| `PullBack(f, D)` |                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `pull_back(f, D)`                                                                                          |
| Takes            | `f: AbstractExpression`, `D: DomainMap \| SideMap \| CellMap` |
| Properties       | `expression` $f$, and `map` $D$, a `DomainMap`, a [`SideMap`](topology.md#sidemap) or a [`CellMap`](topology.md#cellmap) |
| Requires         | $f$ with the domain $Q$ or no domain, and an equal extent for an `Index` free in both $f$ and $D$          |
| Static type      | The value shape of $f$, and the free indices of $f$ followed by those of $D$, with each `Index` taken once |
| Raises           | `ValueError` if another requirement fails                                         |

For a `DomainMap`, the pull-back is the composition $f \circ \psi$, the
expression on $P$ with the value $f(\psi(p))$ at each point $p$ of $P$.
For a `CellMap` or `SideMap`, the map selects cells and reference positions;
its physical coordinates alone do not define the evaluation. The axes of the
value shape are not transformed. An axis of an expression on $Q$ is an array in the
coordinate frame of $Q$, so the pull-back keeps the frame of $Q$.

[`PushForward`](#pushforward) gives the corresponding forward value
transformation from $T$ to the mesh $D$. With the identity law it reverses
the pull-back cellwise; other laws also transform the values. It does not
denote a dual or adjoint operator.

Pull-backs of terminals remain for their supplying libraries to evaluate
(see [Functions in a space](forms.md#functions-in-a-space)).

For coordinate maps $D$ and $E$, with $E$ mapping into the parameter domain
of $D$, their composition is the domain map with the parameterisation
`pull_back(D.parameterisation, E)` on the parameter domain of $E$, in the
ambient domain of $D$, such as `E` in the example below. Composition through
a cell map instead keeps nested pull-backs, preserving the cell selection.

For an edge of the reference triangle `T` and a coordinate map `psi` from
`T` into `RD(2)`:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
one = Integer(1)
e = domain_map(S, tensor((one - t, t)), T)
f = X[0] * X[1]
fe = pull_back(f, e)  # Value (1 - t) * t.
E = domain_map(S, pull_back(psi, e), RD(2))
```

Changing the target to `RD(2)` makes the pull-back of `f` invalid:

```python
edge_in_plane = domain_map(S, tensor((one - t, t)), RD(2))
pull_back(f, edge_in_plane)  # Raises ValueError: f has domain T.
```

### Jacobian and density

The free functions `jacobian(D)`, `inverse_jacobian(D)` and `density(D)` are
the Jacobian matrix, the left inverse of the Jacobian matrix and the density
of a `DomainMap` or a `CellMap` of a mesh, with coordinate parameter domain
$P$ of dimension $k$ and parameterisation $\psi$. $n$ is the physical
coordinate dimension of the target. A `CellMap` of a topology has no
parameterisation and is rejected with `TypeError`.

Each helper checks these input requirements before building its composition.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `jacobian(D)` | | `D: DomainMap \| CellMap` | `gradient(D.parameterisation, P)`, where `P` is `D.parameter_domain` for a `DomainMap` or `D.reference_cell` for a `CellMap` |
| `inverse_jacobian(D)` | | `D: DomainMap \| CellMap` | `inverse(jacobian(D))` |
| `density(D)` | | `D: DomainMap \| CellMap` | The density $\rho$ below |

The Jacobian matrix $J$ has value shape `(n, k)`, and $K$ =
`inverse_jacobian(D)`, of value shape `(k, n)`, is the inverse of $J$ for
$k = n$ and the pseudo-inverse $(J^T J)^{-1} J^T$ for $k < n$. The density
$\rho$ is the factor by which $\psi$ scales the length, the area or the volume
of a small region of $P$:

| Case | `density(D)` |
| ---- | ------------ |
| $k = 0$ | `Integer(1)` |
| $0 < k = n$ | `absolute(determinant(J))` |
| $0 < k < n$ | `square_root(determinant(G))`, with $G = J^T J$ |

The row $k = 0$ gives the value of the other rows, and keeps the Jacobian, of
value shape `(n, 0)`, out of an integral over points.

$J^T J$ is `component_tensor(index_sum(J[l, a] * J[l, b], (l,)), (a, b))`, with
derived indices.

An inverse Jacobian requires $\operatorname{rank} J = k$ wherever evaluated.
The density alone does not require full rank and may vanish.

The core, passes and libraries use `jacobian(D)`, `determinant(jacobian(D))`
and `inverse_jacobian(D)` for $J$, $\det J$ and $K$. These expressions let a
pass identify and evaluate each geometric quantity once per $D$.

### Concept: integral over a domain

Integration of an expression $f$ on $Q$ over $D$ = `DomainMap(P, psi, Q)`
is defined by the parameterisation:

$$
\int_D f = \int_P (f \circ \psi) \, \rho \, \mathrm{d}X,
$$

which is `integral(pull_back(f, D) * density(D), P)` for a scalar $f$. The
integral does not depend on the direction in which $\psi$ traverses $D$.

Coverage is counted with multiplicity. The result equals the ordinary integral
over the geometric image when the parameterisation is injective almost everywhere.

The core has no orientation of a domain. An integral that depends on a
direction, such as a circulation or a flux, is written with the Jacobian
matrix, which carries the direction of $\psi$. For example, the circulation
of a vector field $F$ along an edge $e$ = `DomainMap(S, psi, Q)` with $k = 1$
is

$$
\int_e F \cdot \mathrm{d}r = \int_S (F \circ \psi) \cdot \psi' \, \mathrm{d}t,
$$

which is `integral(index_sum(pull_back(F, e)[i] * jacobian(e)[i, 0], (i,)), S)`,
and the circulation changes sign when $\psi$ runs the other way. Where the
cells of a mesh must be oriented consistently, such as for the
flux through a surface in $R^3$, a separately supplied scalar field gives
a sign $+1$ or $-1$ for each cell (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

### Mapped integration families

An indexed `DomainMap` describes a family of maps and retains its indices.
For a scalar field $f$ on $Q$ and `M = domain_map(P, psi[i, k], Q)`, with a
coordinate parameter domain $P$:

```python
per_member = integral(pull_back(f, M) * density(M), P)
per_configuration = index_sum(per_member, (i,))
```

The first expression retains `i` and `k`. The second sums only over `i` and
retains `k`. For several unrelated maps, add their reference-domain integrals
explicitly. No field domain is created by this construction. Integration counts
each map member separately, including coincident members.

For example, three segment maps of lengths 2, 3 and 5 give an indexed integral
of `Integer(1)` with values 2, 3 and 5. An explicit sum over their member index
gives 10. A topology-backed mesh of those segments also has total length 10,
but additionally identifies a domain on which fields and spaces can be defined.

### Points and point evaluation

A point is the image of a domain map. A point in a domain $Q$ is
`domain_map(RD(0), x, Q)`, with $x$ a constant expression of value shape
`(n,)` for the dimension $n$ of the ambient coordinate domain of $Q$. The core
needs no node class and no free function for a point.

#### Point evaluation

The value of an expression `f` at a point or indexed family of points `p` is:

```python
integral(pull_back(f, p), RD(0))
```

`p` has parameter domain `RD(0)`, and `f` is defined at its points. The
result has no domain, preserves the value shape, scalar type and argument
spaces of `f`, and carries the free indices of `f` followed by those of `p`,
each taken once. Point indices are retained, not summed.

Evaluation preserves symbolic dependence on the expression and the points.
Derivatives and value transformations are applied to the expression before
sampling; sampling `gradient(f, Q)` gives the gradient on `Q` at the points.
Basis tabulation and point observations use the same evaluation. The sampled
values alone define no field between the points.

For points in a mesh, see [A family of points in a
mesh](#a-family-of-points-in-a-mesh).

#### Quadrature

A quadrature rule on a coordinate domain is a weighted sum of point
evaluations. For a family `p` with index `i` and weights `w`:

```python
index_sum(w[i] * integral(pull_back(f, p), RD(0)), (i,))
```

`integral(f, Q)` denotes the integral; this sum approximates it and is
mathematically equal when the rule is exact for `f`. A quadrature pass
applies a supplied rule. On a mapped domain, the integrand includes the density
from [integral expansion](#integral-over-a-mesh).

#### Examples

The examples use the unit triangle `T = UnitSimplex(2)`, `X = coordinate(T)`,
and `f = X[0] * X[1]`:

```python
third = Integer(1) / Integer(3)

# DomainMap: the centroid of T, as a point in T.
c = domain_map(RD(0), tensor((third, third)), T)

# The value at the centroid: 1/9, with no domain and value shape ().
fc = integral(pull_back(f, c), RD(0))
```

Worked examples for this chapter: [an affinely mapped triangle](../examples/integration.md#an-affinely-mapped-triangle), [circulation along an edge](../examples/integration.md#circulation-along-an-edge), [integrating over an edge](../examples/integration.md#integrating-over-an-edge), [integrating over the boundary of a triangle](../examples/integration.md#integrating-over-the-boundary-of-a-triangle), [quadrature](../examples/integration.md#quadrature), [integrating over a cut cell](../examples/integration.md#integrating-over-a-cut-cell), [the centroid of a triangle](../examples/integration.md#the-centroid-of-a-triangle).

## Fields and integration on a mesh

### `TransformationLaw`

`TransformationLaw` specifies how function values on a reference cell transform
to a cell of a mesh.

| `TransformationLaw(variances, weight)` |                                                                       |
| -------------------------------------- | --------------------------------------------------------------------- |
| Extends                                | `Base` |
| Free function                          | `transformation_law(variances, weight)`                               |
| Takes                                  | `variances: tuple[Literal["covariant", "contravariant"], ...]`, `weight: int` |
| Properties                             | `variances`, a `tuple` of `"covariant"` and `"contravariant"`, and `weight` $w$ |

There is one variance for each of the last axes of the value shape, and the
leading axes are unchanged. The permitted variance names are checked statically.

For a cell $\psi$ of a mesh, with $J$, $K$ and $\rho$ of [Jacobian and
density](#jacobian-and-density), $k$ the dimension of the topology and $n$ the
dimension of the ambient coordinate domain, the law $m$ maps a value $\hat f$
on the reference cell to $m_\psi(\hat f)$. $m_\psi$ multiplies each covariant
axis by $K^T$ and each contravariant axis by $J$, and divides by
$(\det J)^w$ for $k = n$ and by $\rho^w$ for $k < n$. Each transformed axis
has the extent $k$ in $\hat f$ and $n$ in $m_\psi(\hat f)$.

The pull-back of a gradient by a cell and the normal of a side use the
covariant law of weight 0 (see [Expanding pull-backs](evaluation.md#expanding-pull-backs)).
The density of an integral is not an instance, because a measure uses
$\lvert \det J \rvert$ and the law uses $\det J$.

The maps of Basix are instances. The core names no instances:

| Space                             | Law                                       | Legacy UFL map         |
| --------------------------------- | ----------------------------------------- | ---------------------- |
| Lagrange, vector Lagrange         | `((), 0)`                                 | identity               |
| Nédélec, $H(\mathrm{curl})$       | `(("covariant",), 0)`                     | covariant Piola        |
| Raviart–Thomas, $H(\mathrm{div})$ | `(("contravariant",), 1)`                 | contravariant Piola    |
| Discontinuous densities           | `((), 1)`                                 | $L^2$ Piola            |
| Regge                             | `(("covariant", "covariant"), 0)`         | double covariant       |
| Hellan–Herrmann–Johnson           | `(("contravariant", "contravariant"), 2)` | double contravariant   |

The law has no orientation. A space of odd weight on a topology of lower
dimension than the ambient coordinate domain gives its own signs for each pair
of cell and facet, as the space gives the transformations of its degrees of
freedom.

### `PushForward`

`PushForward` gives the field on a mesh $D$ of an expression $f$ on its topology
$T$.

| `PushForward(f, D, m)` |                                                                     |
| ---------------------- | ------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Takes                  | `f: AbstractExpression`, `D: Mesh`, `m: TransformationLaw` |
| Properties             | `expression` $f$, `mesh` $D$ and `law` $m$                        |
| Requires               | $D$ a mesh with the topology $T$, $f$ with the domain $T$ or no domain, and the last axes of $f$, one for each variance of $m$, of the extent $k$ |
| Static type            | The value shape of $f$ with the extent $n$ for each transformed axis, the domain $D$ |
| Raises                 | `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `push_forward(f, D, m)` | | `f: AbstractExpression`, `D: Mesh`, `m: TransformationLaw` | `PushForward(f, D, m)` |

The pull-back of `PushForward(f, D, m)` by `cell_map(D, j)` is
$m_\psi(\hat f)$, with $\psi$ the cell and $\hat f$ =
`pull_back(f, cell_map(T, j))`. The push-forward of the coordinate function
$x$ of $D$ by `transformation_law((), 0)` is `coordinate(D)`.

The transformation law is always explicit, for example
`push_forward(u, D, transformation_law(V))` for `u` in `V`, where the
one-argument `transformation_law` is supplied by the library of `V`.
The arguments of a form, the variable and the direction of a derivative,
and the keys of the data are the
functions on $T$. A library may wrap `push_forward`, as the
[classic FEniCS layer](../examples/forms.md#a-classic-fenics-layer-for-the-poisson-problem) does.

The identity of Piola holds on each cell of $D$, with no derivative of $J$:

- The divergence of a push-forward by `(("contravariant",), 1)` is the
  divergence on $K_j$ of $\hat f$, carried by `((), 1)`.
- For $k = n$, the curl of a push-forward by `(("covariant",), 0)` is the curl
  on $K_j$ of $\hat f$, carried by `(("contravariant",), 1)` for $n = 3$ and
  by `((), 1)` for $n = 2$.

[`expand_pull_backs`](evaluation.md#expanding-pull-backs) applies the identity, and every
other `Divergence` and `Curl` lowers to its composition of `Gradient`.

### Traces on a mesh

[Restriction](topology.md#restriction) and [`SideMap`](topology.md#sidemap) select a side on the
topology. On a physical mesh, a side domain also supplies the geometry.

#### Side domains

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_domain(D, S, s)` | | `D: Mesh`, `S: AbstractTopology`, `s: str` | `mesh(pull_back(x, side_map(T, S, s)), Q)` |

For a mesh $D$ = `mesh(x, Q)` with the topology $T$, the side domain of $S$
is the mesh defined by the trace of $x$ from side $s$, so the side domain moves with
$x$. The result has topology $S$ and geometry supplied by that side.

The chosen side supplies the integration measure. Any measure conversion
belongs in the integrand.

#### Traces and normals

The restriction of an expression is the value of the expression on the cell
of $S$, from the side:

- The restriction of `gradient(f, D)` is the whole gradient of `f` on the cell
  on the side, as for a normal derivative. The last axis is in the frame of
  the ambient coordinate domain from every side.
- `gradient(pull_back(f, side_map(D, E, s)), E)` is the tangential gradient of
  the restriction along the cells of $E$.

The unit side conormal is [`normal(side_map(D, E, s))`](#normal), outward on
boundary pieces.

#### Example

The example uses `D` and `square` of the
[`Mesh`](#mesh) example, `I = ReferenceSimplex(1)`, a
function `u` on `square`, and the diagonal of `square`:

```python
diagonal = StaticTopology((I,), ((1, 2),), ((("+", square, 0), ("-", square, 1)),))
plus = side_map(square, diagonal, "+")
minus = side_map(square, diagonal, "-")
jump = pull_back(u, plus) - pull_back(u, minus)
E = side_domain(D, diagonal, "+")
```

The `"+"` side selects cell 0 with inclusion $A_+(t) = (1-t,t)$.
For a reference-cell field $u_0(X) = a_{0,0} + a_{0,1}X_0 + a_{0,2}X_1$,
its trace is $u_+(t) = a_{0,0} + a_{0,1}(1-t) + a_{0,2}t$.
See [Evaluating a trace](../examples/coupling.md#evaluating-a-trace) for both
sides and the effect of reversing the facet vertices.

```{dropdown} Legacy UFL
The legacy UFL restriction `f('+')` is `pull_back(f, side_map(T, S, "+"))` on the topology $S$ of the interior facets, whose sides are named `'+'` and `'-'`.
```

### A family of points in a mesh

The package supplies reference positions in a known local cell slot `j`.
For a domain-free array `positions` of reference positions in $K_j$,
evaluation is:

```python
R = domain_map(RD(0), positions[i], K_j)
values = pull_back(pull_back(f, cell_map(D, j)), R)
```

The point index `i` remains free. It indexes samples within one local cell;
it does not select mesh cells. Point location and partitioning samples among
concrete cells belong to the package. Distinct local slots are read separately.

#### Adjacent and coincident cells

Let a configuration contain two interval slots with reference cell $K = [0, 1]$.
On a mesh $D$ with cells $[0, 1]$ and $[1, 2]$, let $f$ equal 2 on cell 0
and 7 on cell 1. These snippets describe symbolic values, not numerical
execution by the core:

```python
f0 = pull_back(f, cell_map(D, 0))
f1 = pull_back(f, cell_map(D, 1))
midpoint = domain_map(RD(0), tensor((Float(0.5),)), K)
left_endpoint = domain_map(RD(0), tensor((Integer(0),)), K)
right_endpoint = domain_map(RD(0), tensor((Integer(1),)), K)

pull_back(f0, midpoint)  # Physical coordinate 0.5; value 2.
pull_back(f0, right_endpoint)  # Physical coordinate 1; value 2.
pull_back(f1, left_endpoint)  # Physical coordinate 1; value 7.
```

The last two expressions are the two traces at the interface. Physical
coordinate 1 alone does not choose between them.

For a mesh `D_overlap` whose two cells both occupy $[0, 1]$, let $g$ again
have values 2 and 7 on the respective cells:

```python
g0 = pull_back(g, cell_map(D_overlap, 0))
g1 = pull_back(g, cell_map(D_overlap, 1))
pull_back(g0, midpoint)  # Physical coordinate 0.5; value 2.
pull_back(g1, midpoint)  # Physical coordinate 0.5; value 7.
```

The physical and reference coordinates coincide, but the selected cells differ.

#### Points forming a domain

A fixed local collection of points may have a topology `S` with one
zero-dimensional slot per point and a side `("inside", T, j)` selecting
its parent slot. Each slot's side data supply its reference position. Set
`P = side_domain(D, S, "inside")`:

```python
inside = side_map(D, P, "inside")
values = pull_back(pull_back(f, inside), cell_map(P, 0))
total = integral(pull_back(f, inside), P)
```

This supports fields and integration on a fixed point configuration. A
runtime-sized sampling or quadrature array instead uses the indexed
pull-backs above.

### `Integral` over a mesh

`Integral(f, D)` integrates over the fixed local cells of `D`, counting
coincident slots separately. It returns a local contribution. The package
iterates over concrete cells and assembles contributions outside UFLx.
The integrand has domain `D`, domain `D.ambient_domain`, or no domain.
Every free index of `f` remains free.

#### Expansion over cells

For slot $j$ of `D.topology`, let $K_j$ be its reference cell and
`D_j = cell_map(D, j)`. For a scalar field on `D`, its contribution is:

```python
integral(pull_back(f, D_j) * density(D_j), K_j)
```

Expansion adds these terms over the statically known slots. No cell index,
lookup table or cell-axis `IndexSum` is introduced. A single-cell kernel
has one term. A fixed macro-element kernel has one term per subcell.

For a field on `D.ambient_domain`, use the map with the same parameterisation
and that ambient domain as target. Array-valued integrands are integrated
component by component. Zero-dimensional cells have density `Integer(1)`.

The density does not depend on orientation. For ambient fields, the result
agrees with integration over the geometric union when each cell map is injective
almost everywhere and distinct images overlap only on sets of zero measure.
An empty runtime mesh selection requires no kernel calls, not an empty cell sum
inside a kernel.

#### Local data and indexed loads

Local cell-dependent data are fields on the topology or mesh. A field equal
to 2 on one unit cell slot and 7 on another has local integral 9. A domain-free
array with entries `[2, 7]` instead has the same value on both slots; its
integral over that two-unit-cell configuration is `[4, 14]`. Array indices
do not implicitly select slots.

#### Examples

For the local square configuration above, `integral(Integer(1), D)` is the
area 1, the sum of two reference-cell integrals. A boundary configuration
`B` with the side `"inside"` gives the boundary contribution
`integral(pull_back(f, side_map(D, B, "inside")), B)`.

```{dropdown} Legacy UFL
The cell, boundary-facet and interior-facet contributions to `dx`, `ds` and
`dS` are local integrals over their respective configurations. The finite
element package selects concrete entities and assembles the contributions.
```

### `Normal`

`Normal` is the unit side conormal: tangent to the parent cell and orthogonal
to the selected side image in that parent, with orientation supplied by
`reference_normal`.
It points outward on boundary pieces.

| `Normal(sm)` |                                                                                                       |
| ------------ | ----------------------------------------------------------------------------------------------------- |
| Extends      | `AbstractExpression` |
| Free function | `normal(sm)`                                                                                         |
| Takes        | `sm: SideMap` |
| Properties   | `side_map`, a [`SideMap`](topology.md#sidemap) `SideMap(D, E, s)`                                               |
| Requires     | $D$ a mesh, and $E$ of the dimension $k - 1$ for the dimension $k$ of $D$                             |
| Static type  | Value shape `(n,)` for the dimension $n$ of the ambient coordinate domain, the domain $E$, real      |
| Raises       | `ValueError` if the requirement fails           |

The normal of `sm` belongs to the same side as the restriction
`pull_back(f, sm)`, so a form writes both from one side map.

Notation, for `sm` = `SideMap(D, E, s)` and the piece $i$ of $E$:

- $B_i$ and $A_i$ are the maps of [Restriction](topology.md#restriction) for the side $s$.
- $K_j$ is the reference cell of the cell on the side.
- $k$ is the dimension of $D$, and $n$ is the dimension of the ambient
  coordinate domain.

#### The pull-back of `Normal`

The pull-back of `Normal(SideMap(D, E, s))` by `cell_map(E, i)` is
$\nu = w / \lvert w \rvert$, with $w = m_{B_i}(\hat m)$,
using the local maps of [Restriction](topology.md#restriction), where:

- $m$ is `transformation_law(("covariant",), 0)`, applied with the Jacobian
  of `B_i` on the side at the points of `A_i` (see
  [`TransformationLaw`](#transformationlaw)).
- $\hat m$ is the `reference_normal` from `E.topology.side_data(i, s)`
  (see [Side data](topology.md#side-data)).

The composition is `w / square_root(w @ w)`. It requires full column rank
of the cell Jacobian and $w \ne 0$ wherever evaluated.

The core needs no numbering of the facets and no vertices of $K_j$, because
the library of the topology gives $\hat m$.

#### The direction of the normal

Orientation follows the [side-data contract](topology.md#side-data). For $k=n$,
normals are opposite when the two regions adjoin on opposite sides of the
interface; coincident geometry or side names alone do not imply this.
For $k<n$, matching tangent spaces are also required. On a surface, a jump
therefore uses each side's conormal, as in `u_p * n_p + u_m * n_m`.

`Normal` is a side conormal. A unit normal field to a surface in $R^3$
requires a separately supplied orientation field (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

#### Examples

The examples use the following objects:

- `B` of the [`Integral` over a mesh](#integral-over-a-mesh) examples.
- The mesh `E` of the interior facets of `D`, with the sides
  `"+"` and `"-"`.
- `grad_f = gradient(f, D)` and `i = index("i")`.

```python
inside, plus, minus = side_map(D, B, "inside"), side_map(D, E, "+"), side_map(D, E, "-")
n_B = normal(inside)

# The flux of grad f out of the square, the legacy UFL dot(grad(f), n) * ds.
flux = integral(index_sum(pull_back(grad_f, inside)[i] * n_B[i], (i,)), B)

# The jump of the normal derivative across the diagonal, from the side "+".
n_E = normal(plus)
jump = pull_back(grad_f, plus)[i] - pull_back(grad_f, minus)[i]
penalty = integral(index_sum(jump * n_E[i], (i,)) ** Integer(2), E)
```

Worked examples for this chapter: [an interior facet integral](../examples/coupling.md#an-interior-facet-integral), [a Bank-Weiser error estimator](../examples/estimation.md#a-bank-weiser-error-estimator), [a fluid and a solid](../examples/coupling.md#a-fluid-and-a-solid), [a 1D network in a 3D mesh](../examples/coupling.md#a-1d-network-in-a-3d-mesh), [a deformed domain](../examples/surfaces.md#a-deformed-domain).
