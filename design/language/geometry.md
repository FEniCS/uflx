# Geometry and calculus

A local cell configuration acquires physical geometry through a coordinate expression:

```python
# T describes one triangle slot; the package supplies its local data.
Q = EuclideanSpace(2)  # Ambient space: the plane.
X = Lagrange(T, 1, (2,))  # Library-supplied space of coordinate functions.
x = function(X, "x")  # Physical positions, with values supplied by the solver.
D = mesh(x, Q)  # Physical mesh; its topology comes from x.
```

Changing the coordinate expression can change the physical geometry while
retaining the topology and the functions defined on it. Fields on the topology
are carried to the physical mesh by [push-forward](#pushforward).

The only metric is the ambient Euclidean metric, which defines the tangential
[gradient](#gradient), integration [density](#jacobian-and-density) and unit
[normal](#normal).

## Meshes

### `Mesh`

`Mesh(coordinate_function, ambient_domain)` gives a local cell configuration
its physical geometry. `coordinate_function` supplies physical coordinates
on the topology; `ambient_domain` is the coordinate domain containing them.
Cells keep their identities when their physical coordinates coincide. The mesh stores no
global connectivity or coordinate values. The package supplies the local data.
Its finitely many cell slots describe one kernel configuration, not a runtime mesh.

| `Mesh(coordinate_function, ambient_domain)` | |
| ------------ | - |
| Extends | `AbstractDomain` |
| Free function | `mesh(coordinate_function, ambient_domain)` |
| Takes | `coordinate_function: AbstractExpression`, `ambient_domain: AbstractCoordinateDomain` |
| Properties | `coordinate_function`, the coordinate expression on the topology, and `ambient_domain`, the coordinate domain containing its physical points |
| Requires | `coordinate_function` has a topology domain whose cell types are reference cells; it is real, has value shape `(n,)` with `n = ambient_domain.dimension`, and has free indices `()` and argument spaces `()`; `topology.dimension <= n` |
| Derived | `topology` is `coordinate_function.domain`; `dimension` is `topology.dimension` |
| Raises | `TypeError` if the coordinate expression has no topology domain, and `ValueError` if another requirement fails |

The factory constructs the corresponding `Mesh` directly. The topology is
determined by `coordinate_function`; it is not a separate construction
argument. Mesh identity depends on both stored properties. Replacing the
coordinates reconstructs the mesh from the new coordinates.

In the descriptions below, `D = mesh(x, Q)`, `x` is the coordinate expression,
`Q` is its ambient domain, and `T = x.domain` is its topology.

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

#### Mesh cell maps

For a mesh `D = mesh(x, Q)` on topology `T`, `cell_map(D, j)` selects local
slot `j`. Its reference cell is `T.cell_types[j]` and its derived
parameterisation is `pull_back(x, cell_map(T, j))`. Pulling back
`coordinate(D)` gives this same expression.

The cell map retains the selected cell when its coordinates change:
substitution reconstructs it from `D` and `j`. For a field `f` on `D`,
`pull_back(f, cell_map(D, j))` is its reference-cell representation.

To read a field on the ambient coordinate domain instead, construct a
coordinate map from the cell map:

```python
cell = cell_map(D, j)
coordinate_map = domain_map(cell.reference_cell, cell.parameterisation, D.ambient_domain)
```

These maps have the same parameterisation but different targets: `cell`
selects a cell of `D`, while `coordinate_map` gives points in `D.ambient_domain`.
A mesh embedded in another mesh still has its own physical coordinates;
cell or side maps express the relation to the other mesh's cells.

#### Examples

Use the two-triangle topology of the [static topology
example](../examples/library.md#a-static-topology), on reference cells:

```python
K = ReferenceSimplex(2)
square = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)))
```

Let `x` be a vector Lagrange function of degree 1 on `square`, of value shape
`(2,)`, with vertex values $(0, 0)$, $(1, 0)$, $(0, 1)$ and $(1, 1)$.
The resulting physical mesh is the square:

```python
D = mesh(x, EuclideanSpace(2))
```

- `D.topology` is `square`, with two triangle slots. `cell_map(D, 0)`
  reads the first; `cell_map(D, 1)` reads the second.
- `pull_back(coordinate(D), cell_map(D, 1))` is
  `pull_back(x, cell_map(square, 1))`, which the library evaluates to
  $(1 - X_1, X_0 + X_1)$ for the supplied coordinates.

For a local configuration `M` with types `(Kt, Kq)`, `mesh(x, EuclideanSpace(2))` has
one triangle cell and one quadrilateral cell. A field `u` on `M` has local reads:

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

| `Coordinate(domain)` |                                                    |
| --------------- | -------------------------------------------------- |
| Extends         | `AbstractExpression` |
| Free function   | `coordinate(domain)`                                    |
| Takes           | `domain: AbstractDomain` |
| Properties | `domain`, the domain whose point coordinates are read |
| Static type | Value shape `(n,)`, real, with domain `domain`; `n` is the dimension of its ambient coordinate domain |

Let $P$ denote `domain` and $n$ its ambient coordinate dimension. On a
coordinate domain, the coordinate function is the identity,
$X : P \to R^n$, $X(p) = p$, and the gradient of $X$ is the identity matrix
(see [`Gradient`](#gradient)).

For a mesh $D$, `coordinate(D)` is the physical coordinate field on $D$.
Its pull-back by a cell map is that cell's parameterisation (see
[`Mesh`](#mesh)). This differs from `D.coordinate_function`, which is
the stored coordinate expression on the topology.
For a mesh `mesh(x, Q)`, the pull-back by `cell_map(D, j)` is
`pull_back(x, cell_map(T, j))`.

For example, `coordinate(EuclideanSpace(3))` has value shape `(3,)` and domain `EuclideanSpace(3)`.

### `Gradient`

`Gradient(expression, domain)` differentiates `expression` with respect to
position on `domain`. It adds a last array axis, with one entry for each
ambient coordinate direction. Let `n` be the dimension of the ambient
coordinate domain.

| `Gradient(expression, domain)` |                                                                                                   |
| ---------------- | ------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `gradient(expression, domain)`                                                                                  |
| Takes            | `expression: AbstractExpression`, `domain: AbstractDomain` |
| Properties | `expression`, the field to differentiate, and `domain`, where the spatial derivative is taken |
| Requires         | `expression` has domain `domain` or no domain                                                              |
| Static type      | Value shape `expression.value_shape + (n,)`                                                                |
| Raises           | `ValueError` for `expression` if a requirement fails                               |

Writing $f$ for the expression and $P$ for its domain, the gradient acts on
each component of an array-valued $f$. Entry `j` of the
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

Start with a scalar expression on `EuclideanSpace(3)`:

```python
P = EuclideanSpace(3)
X = coordinate(P)
x, y, z = X[0], X[1], X[2]
f = x * y * z
```

| Expression                              | Value shape | Meaning                                       |
| --------------------------------------- | ----------- | --------------------------------------------- |
| `gradient(f, P)`                        | `(3,)`      | The gradient $(yz, xz, xy)$, after evaluation |
| `gradient(gradient(f, P), P)`           | `(3, 3)`    | The Hessian of `f`                            |
| `gradient(X, P)`                        | `(3, 3)`    | The identity matrix, after evaluation         |

- Let `M` be a mesh of dimension 2 in `EuclideanSpace(3)`, and `w` a field
  on `M`. `gradient(w, M)` is the surface gradient of `w`, with value shape
  `(3,)`. `gradient(coordinate(M), M)` is the orthogonal projection onto the
  tangent plane of `M`, with value shape `(3, 3)`.
- `gradient(X, EuclideanSpace(2))` raises `ValueError`, because `X` has the domain
  `EuclideanSpace(3)`.

### `Divergence` and `Curl`

`Divergence(expression, domain)` contracts the last value axis with the
spatial derivative axis. `Curl(expression, domain)` takes the curl of a vector
field in two or three ambient dimensions. Let `n` be the dimension of the
ambient coordinate domain.

| `Divergence(expression, domain)` |                                                                     |
| ------------------ | ------------------------------------------------------------------- |
| Extends            | `AbstractExpression` |
| Free function      | `divergence(expression, domain)`                                                  |
| Takes              | `expression: AbstractExpression`, `domain: AbstractDomain` |
| Properties | `expression`, the field to differentiate, and `domain`, where the spatial derivative is taken |
| Requires           | `expression` has domain `domain` or no domain, and `expression` of rank at least 1 with last extent `n` |
| Static type        | Value shape `expression.value_shape[:-1]`                                    |
| Raises             | `ValueError` for `expression` if a requirement fails |

| `Curl(expression, domain)`  |                                                                     |
| ------------- | ------------------------------------------------------------------- |
| Extends       | `AbstractExpression` |
| Free function | `curl(expression, domain)`                                                        |
| Takes         | `expression: AbstractExpression`, `domain: AbstractDomain` |
| Properties | `expression`, the field to differentiate, and `domain`, where the spatial derivative is taken |
| Requires      | `expression` has domain `domain` or no domain, `n` is 2 or 3, and `expression` of value shape `(n,)` |
| Static type   | Value shape `(3,)` for `n == 3` and `()` for `n == 2`                 |
| Raises        | `ValueError` for `expression` if a requirement fails |

Writing `u` for `expression` and `P` for `domain`, the value of
`Divergence(u, P)` is the trace of `gradient(u, P)` over the last
axis of `u` and the axis of the gradient. With `G = gradient(u, P)`, so
$G_{ij}=\partial_j u_i$, `Curl(u, P)` is the antisymmetric contraction:

| Ambient dimension | Value |
| ----------------- | ----- |
| 2 | `G[1, 0] - G[0, 1]` |
| 3 | `tensor((G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]))` |

`Gradient` stays the one derivative that defines the others. The two node
classes are in the core so that the identity of Piola
holds structurally (see [`PushForward`](#pushforward)).

### `DirectionalDerivative`

`DirectionalDerivative(expression, variable, point, direction)` varies the
specified function variable and measures the first-order change in
`expression`. `point` supplies the variable's value at which to differentiate;
`direction` supplies its perturbation. Here a point is a value in the
variable's function space, not a spatial location. This is the *Gâteaux
directional derivative*, taken with a real perturbation parameter. The node
represents its value for the supplied direction.

Writing $f$, $x$, $a$ and $h$ for these four arguments, its value is

$$
\mathrm{D}(f, x, a, h) = \frac{d}{d\varepsilon} f[x := a + \varepsilon h] \Big|_{\varepsilon = 0},
$$

with $\varepsilon$ real. Substitute free occurrences of $x$ in a compatible
fixed-space representation of $f$ with [geometry dependence](#geometry-dependence) exposed.
For tuple variables, matching tuples $a$ and $h$ supply all replacements
$x_k:=a_k+\varepsilon h_k$ simultaneously. Members may mix real and complex
spaces. The functions of $x$ are bound in $f$, not in $a$ or $h$.

| `DirectionalDerivative(expression, variable, point, direction)` |                                                                                                       |
| -------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Extends                          | `AbstractExpression` |
| Takes                            | `expression: AbstractExpression`, `variable: Function \| tuple[Function, ...]`, `point, direction: AbstractExpression \| tuple[AbstractExpression, ...]` |
| Properties | `expression`, the body to differentiate; `variable`, the function or tuple of functions bound in that body; `point`, their evaluation values; and `direction`, their perturbations |
| Requires                         | `variable` is a [bindable function](forms.md#bindable-functions), or a `tuple` of at least two distinct bindable functions; member by member, `point` and `direction` admissible values for the space of `variable` as in [argument admissibility](forms.md#admissible-argument-values); variable spaces independent of all variable members |
| Static type                      | The value shape of `expression`, the free indices of `expression` followed by those of the members of `point` and `direction`, each `Index` taken once, and the argument spaces and the scalar type of `expression` |
| Raises                           | `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `directional_derivative(expression, variable, point, direction)` | | `expression: AbstractExpression`, `variable: Function \| tuple`, `point, direction: AbstractExpression \| tuple` | `DirectionalDerivative(expression, variable, point, direction)` |

A real-linear differential is assumed, even for complex variables and
directions. Real points and directions enter complex spaces by
[complex embedding](rewriting.md#type-compatibility-of-rewrites).
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
See [Geometry preparation](differentiation.md#geometry-preparation).

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

`Integral(expression, integration_domain)` integrates `expression` over
`integration_domain`, component by component. The result has no spatial
domain, even though it may depend on the integrand and the domain's geometry.
The two core cases are coordinate domains and meshes.

| `Integral(expression, integration_domain)` |                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `integral(expression, integration_domain)`                                                                            |
| Takes            | `expression: AbstractExpression`, `integration_domain: AbstractDomain` |
| Properties | `expression`, the integrand, and `integration_domain`, the region over which it is integrated |
| Requires | `expression` has domain `integration_domain` or no domain; a mesh also accepts an expression on `integration_domain.ambient_domain` |
| Static type      | The value shape and free indices of `expression`, with no domain |
| Raises           | `ValueError` for `expression` if a requirement fails                         |

#### Integration over coordinate domains

Let $f$ denote the integrand and $P$ the integration domain. When $P$ is a
coordinate domain of dimension $n$, `Integral(f, P)` is the Lebesgue integral
in the coordinates $X$ of $P$,

$$
\int_P f \,\mathrm{d}X.
$$

For $n = 0$, the integral is the value of $f$ at the single point of $R^0$.

A later pass may approximate a coordinate-domain integral by
[quadrature](#quadrature), or evaluate it exactly when the rule is exact.
Expressions with no domain have the
[constant-field interpretation](expressions.md#domain-free-expressions).

For example, with `T = UnitSimplex(2)` and `X = coordinate(T)`:

- `integral(X[0], T)` is $\int_T X_0 \,\mathrm{d}X$. A later pass evaluates
  the integral to $1/6$.
- `integral(Integer(1), T)` is the area $1/2$ of `T`.
- `integral(X[0], T) * integral(X[1], T)` is a product of two integrals, with
  value shape `()` and no domain.

#### Integration over mesh cells

`Integral(f, D)` integrates over the fixed local cells of `D`, counting
coincident slots separately. It returns a local contribution. The package
iterates over concrete cells and assembles contributions outside UFLx.
For slot $j$ of `D.topology`, let $K_j$ be its reference cell and
`D_j = cell_map(D, j)`. For a scalar field on `D`, its contribution is:

```python
integral(pull_back(f, D_j) * density(D_j), K_j)
```

Expansion adds these terms over the statically known slots. No cell index,
lookup table or cell-axis `IndexSum` is introduced. A single-cell kernel
has one term. A fixed macro-element kernel has one term per subcell.

For a field on `D.ambient_domain`, replace `D_j` in the pull-back by
`domain_map(K_j, D_j.parameterisation, D.ambient_domain)`. Array-valued
integrands are integrated component by component. Zero-dimensional cells have density `Integer(1)`.

The [density](#jacobian-and-density) does not depend on orientation. For ambient fields, the result
agrees with integration over the geometric union when each cell map is injective
almost everywhere and distinct images overlap only on sets of zero measure.
An empty runtime mesh selection requires no kernel calls, not an empty cell sum
inside a kernel.

Local cell-dependent data are fields on the topology or mesh. A field equal
to 2 on one unit cell slot and 7 on another has local integral 9. A domain-free
array with entries `[2, 7]` instead has the same value on both slots; its
integral over that two-unit-cell configuration is `[4, 14]`. Array indices
do not implicitly select slots.

For the local square configuration above, `integral(Integer(1), D)` is the
area 1, the sum of two reference-cell integrals. A boundary configuration
`B` with the side `"inside"` gives the boundary contribution
`integral(pull_back(f, side_map(D, B, *boundary_sides["inside"].entries)), B)`.

```{dropdown} Legacy UFL
The cell, boundary-facet and interior-facet contributions to `dx`, `ds` and
`dS` are local integrals over their respective configurations. The finite
element package selects concrete entities and assembles the contributions.
```

## Maps between domains

### `DomainMap`

`DomainMap(parameter_domain, parameterisation, ambient_domain)` describes a
geometric map by giving the target coordinates at each parameter point.
`parameter_domain` is where the map starts; `parameterisation` computes the
coordinates; `ambient_domain` is the coordinate domain containing its image.

Writing $P$, $\psi$ and $Q$ for these arguments, the map is $\psi: P \to Q$
and its geometric image is $\psi(P)$. The sections below write $D$ for the
map; its image is always named separately. A `DomainMap` is not an
`AbstractDomain` and is not the domain of an expression.

| `DomainMap(parameter_domain, parameterisation, ambient_domain)` |                                                                                                                                                                                                             |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                          | `Base` |
| Free function                    | `domain_map(parameter_domain, parameterisation, ambient_domain)`                                                                                                                                                                           |
| Takes                            | `parameter_domain: AbstractCoordinateDomain \| AbstractTopology`, `parameterisation: AbstractExpression`, `ambient_domain: AbstractCoordinateDomain` |
| Properties | `parameter_domain`, the domain of the map inputs; `parameterisation`, the coordinate expression; and `ambient_domain`, the coordinate domain of the outputs |
| Requires | `parameterisation` has domain `parameter_domain` or no domain; it is real, has value shape `(n,)` with `n = ambient_domain.dimension`, and has argument spaces `()`; `parameter_domain.dimension <= n` |
| Derived                          | `free_indices`, the free indices of `parameterisation` |
| Raises                           | `ValueError` if another requirement fails                                                                                                                                   |

The topological dimension is `parameter_domain.dimension`; the geometric
dimension is `ambient_domain.dimension`. Maps need not be injective.
Coordinates alone do not select a cell or reference preimage in a mesh;
mesh-field evaluation uses [cell or side maps](#points-in-a-mesh).

A parameterisation with free indices makes `D` a family of maps, with one
map for each assignment of those indices. For example, with an array `x` of value shape
`(m, n)` and an index `i`, `domain_map(EuclideanSpace(0), x[i], EuclideanSpace(n))` is a family of $m$
points in `EuclideanSpace(n)`, retaining distinct members even at identical coordinates.

For example, the hypotenuse of the unit triangle has parameter domain
$P = [0,1]$ and parameterisation $\psi(t) = (1-t,t)$. Choose `Q = EuclideanSpace(2)`
for the edge in the plane, or `Q = UnitSimplex(2)` for the edge in the triangle.

### `PullBack`

`PullBack(expression, map)` reads `expression` through `map`. The expression
is defined on the map's target, or has no domain; the result is defined on
the domain from which the map starts. Pull-back changes where a field is
read, without transforming its value axes.

In the formulas below, $f$ denotes the expression, $D$ the map and $\psi$
its coordinate parameterisation when it is a `DomainMap`.

Here $P$ and $Q$ denote the domains from which and into which the map goes.
Their API properties depend on the map's role:

| Map class | Starting domain $P$ | Target $Q$ |
| --------- | ------------------- | ---------- |
| `DomainMap` | `parameter_domain` | `ambient_domain` |
| `CellMap` | `reference_cell` | `parent` |
| `SideMap` | `source` | `parent` |

| `PullBack(expression, map)` |                                                                                                            |
| ---------------- | ---------------------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `pull_back(expression, map)`                                                                                          |
| Takes            | `expression: AbstractExpression`, `map: DomainMap \| SideMap \| CellMap` |
| Properties | `expression`, the field to read, and `map`, the coordinate, cell or side map through which it is read |
| Requires | `expression` has the target domain listed above or no domain; any `Index` free in both `expression` and `map` has equal extent in both |
| Static type      | The value shape of `expression`, and the free indices of `expression` followed by those of `map`, with each `Index` taken once |
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
`T` into `EuclideanSpace(2)`:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
one = Integer(1)
e = domain_map(S, tensor((one - t, t)), T)
f = X[0] * X[1]
fe = pull_back(f, e)  # Value (1 - t) * t.
E = domain_map(S, pull_back(psi, e), EuclideanSpace(2))
```

Changing the target to `EuclideanSpace(2)` makes the pull-back of `f` invalid:

```python
edge_in_plane = domain_map(S, tensor((one - t, t)), EuclideanSpace(2))
pull_back(f, edge_in_plane)  # Raises ValueError: f has domain T.
```

### Jacobian and density

The free functions `jacobian(map)`, `inverse_jacobian(map)` and `density(map)` are
the Jacobian matrix, the left inverse of the Jacobian matrix and the density
of a `DomainMap` or a `CellMap` of a mesh, with coordinate parameter domain
$P$ of dimension $k$ and parameterisation $\psi$. $n$ is the physical
coordinate dimension of the target. A `CellMap` of a topology has no
parameterisation and is rejected with `TypeError`.

Each helper checks these input requirements before building its composition.
The formulas below write `D` for `map`.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `jacobian(map)` | | `map: DomainMap \| CellMap` | `gradient(map.parameterisation, P)`, where `P` is `map.parameter_domain` for a `DomainMap` or `map.reference_cell` for a `CellMap` |
| `inverse_jacobian(map)` | | `map: DomainMap \| CellMap` | `inverse(jacobian(map))` |
| `density(map)` | | `map: DomainMap \| CellMap` | The density $\rho$ below |

The Jacobian matrix $J$ has value shape `(n, k)`, and $K$ =
`inverse_jacobian(D)`, of value shape `(k, n)`, is the inverse of $J$ for
$k = n$ and the pseudo-inverse $(J^T J)^{-1} J^T$ for $k < n$. The density
$\rho$ is the factor by which $\psi$ scales the length, the area or the volume
of a small region of $P$:

| Case | `density(map)` |
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

### Integration through a map

For a scalar field `f` on `Q` and a map `D = domain_map(P, psi, Q)` with a
coordinate parameter domain `P`, integration through the map is the composition

```python
integral(pull_back(f, D) * density(D), P)
```

Writing $\psi$ for the parameterisation and $\rho$ for its density, this is

$$
\int_P (f \circ \psi)\,\rho\,\mathrm{d}X.
$$

`Integral` still takes the domain `P`, not the map `D`. Coverage is counted
with multiplicity; this equals integration over the image $\psi(P)$ when
the parameterisation is injective almost everywhere. Reversing the traversal
does not change the density.

Direction-dependent integrals use the Jacobian explicitly. For a vector
field `F` and an edge map `edge` with parameter domain `S`, circulation is

```python
integral(index_sum(pull_back(F, edge)[i] * jacobian(edge)[i, 0], (i,)), S)
```

Reversing the edge reverses its tangent and hence the circulation. A flux
through an oriented surface likewise needs a separately supplied orientation
field; see [The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface).

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

### Evaluation at points

Point evaluation uses a `DomainMap` whose parameter domain is
`EuclideanSpace(0)`. For a point or indexed family of points represented by
`point_map`, the common operation is

```python
integral(pull_back(field, point_map), EuclideanSpace(0))
```

The result has no spatial domain. It preserves the field's value shape,
scalar type and argument spaces, and carries the free indices of `field`
followed by those of `point_map`, each taken once. Point indices remain free;
only an explicit sum combines samples.

#### Points in a coordinate domain

For a coordinate domain `Q`, use `domain_map(EuclideanSpace(0), position, Q)`.
`position` is a domain-free expression with value shape `(Q.dimension,)`;
free indices describe a family of points. No separate point class is needed.

#### Points in a mesh

A mesh-field value additionally needs a cell selection. The package locates
points in concrete cells and supplies reference positions within a known
local slot. For a mesh `D`, local slot `j` and a domain-free array `positions`
of reference coordinates:

```python
cell = cell_map(D, j)
local_field = pull_back(field, cell)
point_map = domain_map(EuclideanSpace(0), positions[i], cell.reference_cell)
values = integral(pull_back(local_field, point_map), EuclideanSpace(0))
```

`i` selects samples within that cell; it does not select mesh cells. Each
local slot is read separately. For example, adjacent intervals `[0, 1]` and
`[1, 2]` carrying the constant values 2 and 7 have different traces at
physical coordinate 1. The selected cell determines which value is read,
even when cells overlap. See [Adjacent and coincident cells](../examples/integration.md#adjacent-and-coincident-cells).

A fixed point configuration that must support fields and integration can
instead be a mesh of zero-dimensional cells, with side-map entries selecting
its parents. See [Points forming a domain](../examples/integration.md#points-forming-a-domain).

Evaluation preserves symbolic dependence on both fields and points. Apply
spatial derivatives and value transformations before sampling; for example,
sample `gradient(field, Q)` to obtain its gradient on `Q`. Basis tabulation
and point observations use the same operation. Samples alone define no field
between points.

#### Quadrature

A quadrature rule on a coordinate domain is a weighted sum of point
evaluations. For a point-map family `p` with index `i` and weights `w`:

```python
index_sum(w[i] * integral(pull_back(f, p), EuclideanSpace(0)), (i,))
```

This approximates `integral(f, Q)` and is mathematically equal when the rule
is exact for `f`. A quadrature pass applies a supplied rule. Integration
through a map includes its density in the sampled integrand.

Worked examples: [quadrature](../examples/integration.md#quadrature),
[cut-cell integration](../examples/integration.md#integrating-over-a-cut-cell),
and [point evaluation in a triangle](../examples/integration.md#point-evaluation-in-a-triangle).

## Fields on a mesh

### `TransformationLaw`

`TransformationLaw` specifies how function values on a reference cell transform
to a cell of a mesh.

| `TransformationLaw(variances, weight)` |                                                                       |
| -------------------------------------- | --------------------------------------------------------------------- |
| Extends                                | `Base` |
| Free function                          | `transformation_law(variances, weight)`                               |
| Takes                                  | `variances: tuple[Literal["covariant", "contravariant"], ...]`, `weight: int` |
| Properties                             | `variances`, the transformation for each trailing value axis, and `weight`, the integer exponent of the determinant or density divisor |

There is one variance for each of the last axes of the value shape, and the
leading axes are unchanged. The permitted variance names are checked statically.

Let $w$ denote `weight`. For a cell $\psi$ of a mesh, with $J$, $K$ and
$\rho$ of [Jacobian and density](#jacobian-and-density), $k$ the dimension of the topology and $n$ the
dimension of the ambient coordinate domain, the law $m$ maps a value $\hat f$
on the reference cell to $m_\psi(\hat f)$. $m_\psi$ multiplies each covariant
axis by $K^T$ and each contravariant axis by $J$, and divides by
$(\det J)^w$ for $k = n$ and by $\rho^w$ for $k < n$. Each transformed axis
has the extent $k$ in $\hat f$ and $n$ in $m_\psi(\hat f)$.

The pull-back of a gradient by a cell and the normal of a side use the
covariant law of weight 0 (see [Expanding pull-backs](lowering.md#expanding-pull-backs)).
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

`PushForward(expression, mesh, law)` carries a field from the mesh's topology
to the physical mesh. `law` specifies how its value axes and density weight
transform under each cell's coordinate map.

| `PushForward(expression, mesh, law)` |                                                                     |
| ---------------------- | ------------------------------------------------------------------- |
| Extends                | `AbstractExpression` |
| Takes                  | `expression: AbstractExpression`, `mesh: Mesh`, `law: TransformationLaw` |
| Properties | `expression`, the field on the topology; `mesh`, the physical destination; and `law`, the transformation of field values |
| Requires | `expression` has domain `mesh.topology` or no domain; each trailing value axis selected by `law.variances` has extent `mesh.topology.dimension` |
| Static type | The value shape of `expression` with each transformed axis replaced by the ambient coordinate dimension of `mesh`; domain `mesh` |
| Raises                 | `ValueError` if another requirement fails |

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `push_forward(expression, mesh, law)` | | `expression: AbstractExpression`, `mesh: Mesh`, `law: TransformationLaw` | `PushForward(expression, mesh, law)` |

Writing `f`, `D` and `m` for these arguments and `T = D.topology`, let $k$
be the topology dimension and $n$ the ambient coordinate dimension.
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

[`expand_pull_backs`](lowering.md#expanding-pull-backs) applies the identity, and every
other `Divergence` and `Curl` lowers to its composition of `Gradient`.

### Traces on a mesh

[Restriction](topology.md#restriction) evaluates through an explicit
[`SideMap`](topology.md#sidemap). Its source mesh supplies integration geometry.

#### Side domains

The optional helper constructs source geometry from a chosen relation. This
choice is independent of subsequent restrictions into other parents.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_domain(parent, source, *entries)` | | `parent: Mesh`, `source: AbstractTopology`, `entries: tuple[int, DomainMap]` | `mesh(pull_back(parent.coordinate_function, side_map(parent.topology, source, *entries)), parent.ambient_domain)` |

Coordinate functions may be vector-valued or tensors of scalar functions. For a
conforming interface, the package supplies coordinates and facet correspondence
whose physical traces agree; the coordinate element alone does not ensure this.

#### Traces and normals

The restriction of an expression is the value of the expression on the cell
of $S$, from the side:

- The restriction of `gradient(f, D)` is the whole gradient of `f` on the cell
  on the side, as for a normal derivative. The last axis is in the frame of
  the ambient coordinate domain from every side.
- `gradient(pull_back(f, side_map(D, E, *diagonal_sides[s].entries)), E)` is the tangential gradient of
  the restriction along the cells of $E$.

The unit side conormal is [Normal](#normal), with orientation supplied
independently of the evaluation map.

#### Example

The example uses `D` and `square` of the
[`Mesh`](#mesh) example, `I = ReferenceSimplex(1)`, a
function `u` on `square`, and the diagonal of `square`:

```python
diagonal = StaticTopology((I,), ((1, 2),))
t = coordinate(I)[0]
A_plus = domain_map(I, tensor((Integer(1) - t, t)), K)
A_minus = domain_map(I, tensor((Integer(0), t)), K)
plus = side_map(square, diagonal, (0, A_plus))
minus = side_map(square, diagonal, (1, A_minus))
jump = pull_back(u, plus) - pull_back(u, minus)
E = side_domain(D, diagonal, (0, A_plus))
```

The map `plus` selects cell 0 with inclusion $A_+(t) = (1-t,t)$.
For a reference-cell field $u_0(X) = a_{0,0} + a_{0,1}X_0 + a_{0,2}X_1$,
its trace is $u_+(t) = a_{0,0} + a_{0,1}(1-t) + a_{0,2}t$.
See [Evaluating a trace](../examples/coupling.md#evaluating-a-trace) for both
sides and the effect of reversing the facet vertices.

```{dropdown} Legacy UFL
The legacy UFL restriction `f('+')` is `pull_back(f, side_map(T, S, "+"))` on the topology $S$ of the interior facets, whose sides are named `'+'` and `'-'`.
```

### `Normal`

`Normal` transports and normalises explicitly oriented reference conormals.
Changing orientation leaves the side map and its pull-backs unchanged.

| `Normal(side_map, reference_conormals)` | |
| ------------------------------------ | - |
| Extends | `AbstractExpression` |
| Takes | `side_map: SideMap`, `reference_conormals: tuple[AbstractExpression, ...]` |
| Free function | `normal(side_map, *reference_conormals)` |
| Requires | A mesh parent; codimension one; one real conormal per source slot, without arguments or free indices |
| Static type | Vector of the parent's ambient dimension, on `side_map.source`; no free indices or argument spaces |
| Raises | `ValueError` for an invalid count, shape, domain, scalar type or free indices |

For entry `(j, A_i)`, conormal $\hat m_i$ has shape `(K_j.dimension,)` and
domain `A_i.parameter_domain` or `None`. The library supplies a nonzero vector
orthogonal to the reference inclusion with its chosen orientation. Neither
orthogonality nor nondegeneracy is established by shape validation.

```python
outward = normal(plus, m_plus)
inward = normal(plus, -m_plus)
# Both expressions use the identical evaluation map `plus`.
```

#### The pull-back of Normal

For `sm.entries[i] = (j, A_i)`, let `B_i = cell_map(sm.parent, j)` and
$K = J^+$, the left inverse of the parent cell Jacobian. Pulling back
`Normal(sm, reference_conormals)` by `cell_map(sm.source, i)` gives

$$
\frac{w}{|w|}, \qquad
w = \operatorname{pull\_back}(K^T, A_i)\,\hat m_i.
$$

This requires a full-column-rank cell Jacobian and $w\ne0$. The vector lies
in the parent's tangent space and is orthogonal to the mapped facet. The source
geometry does not choose its orientation.

#### The direction of the normal

Orientation follows the supplied reference conormals. For $k=n$,
normals are opposite when the two regions adjoin on opposite sides of the
interface; coincident geometry or side names alone do not imply this.
For $k<n$, matching tangent spaces are also required. On a surface, a jump
therefore uses each side's conormal, as in `u_p * n_p + u_m * n_m`.

`Normal` is a side conormal. A unit normal field to a surface in $R^3$
requires a separately supplied orientation field (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

#### Examples

The examples use the following objects:

- `B` of the [integration over mesh cells](#integration-over-mesh-cells) examples.
- The mesh `E` of the interior facets of `D`, with the sides
  `"+"` and `"-"`.
- `grad_f = gradient(f, D)` and `i = index("i")`.

```python
inside, plus, minus = side_map(D, B, *boundary_sides["inside"].entries), side_map(D, E, *diagonal_sides["+"].entries), side_map(D, E, *diagonal_sides["-"].entries)
n_B = outward_normal(inside)

# The flux of grad f out of the square, the legacy UFL dot(grad(f), n) * ds.
flux = integral(index_sum(pull_back(grad_f, inside)[i] * n_B[i], (i,)), B)

# The jump of the normal derivative across the diagonal, from the side "+".
n_E = outward_normal(plus)
jump = pull_back(grad_f, plus)[i] - pull_back(grad_f, minus)[i]
penalty = integral(index_sum(jump * n_E[i], (i,)) ** Integer(2), E)
```

Worked examples for this chapter: [an interior facet integral](../examples/coupling.md#an-interior-facet-integral), [a Bank-Weiser error estimator](../examples/estimation.md#a-bank-weiser-error-estimator), [a fluid and a solid](../examples/coupling.md#a-fluid-and-a-solid), [a 1D network in a 3D mesh](../examples/coupling.md#a-1d-network-in-a-3d-mesh), [a deformed domain](../examples/surfaces.md#a-deformed-domain).
