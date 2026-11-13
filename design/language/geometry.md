# Geometry and calculus

## Meshes

### `Mesh`

`Mesh(coordinate_function, ambient_domain)` describes cells whose physical
coordinates are given by `coordinate_function` in `ambient_domain`. The domain
of `coordinate_function` is the topology that specifies the cell slots.

A `Mesh` describes a fixed local cell configuration, such as a cell or an
interface pair, with physical geometry supplied by the coordinate field.

| `Mesh(coordinate_function, ambient_domain)` | |
| ------------ | - |
| Extends | `AbstractDomain` |
| Free function | `mesh(coordinate_function, ambient_domain)` |
| Takes | `coordinate_function: AbstractExpression`, `ambient_domain: AbstractCoordinateDomain` |
| Properties | `coordinate_function`, the coordinate expression on the topology, and `ambient_domain`, the coordinate domain containing its physical points |
| Requires | `coordinate_function` has a topology domain whose cell types are reference cells; it is real, has value shape `(n,)` with `n = ambient_domain.dimension`, and has free indices `()` and argument spaces `()`; `topology.dimension <= n` |
| Derived | `topology` is `coordinate_function.domain`; `dimension` is `topology.dimension` |
| Raises | `TypeError` if the coordinate expression has no topology domain, and `ValueError` if another requirement fails |

Mesh identity depends on both properties. Replacing coordinates reconstructs
the mesh and preserves its requirements. Below, `physical_mesh = mesh(coordinates, ambient_domain)` and
`topology = coordinates.domain`.

#### Mesh cell maps

For `physical_mesh = mesh(coordinates, ambient_domain)`, let
`topology = coordinates.domain`. Indexing `physical_mesh[cell_slot]` selects
that local cell. Its reference cell is `topology.cell_types[cell_slot]`, and
its parameterisation is `pull_back(coordinates, topology[cell_slot])`.

The cell map retains its parent and slot when coordinates change. Lowering
uses `pull_back(u, physical_mesh[cell_slot])` to express a mesh field on the
selected reference cell.

#### Examples

Use the two-triangle topology of the [static topology
example](../examples/library.md#a-static-topology), on reference cells:

```python
reference_triangle = ReferenceSimplex(2)
square_topology = StaticTopology((reference_triangle, reference_triangle), ((0, 1, 2), (1, 3, 2)))
```

Let `coordinates` be a vector Lagrange function of degree 1 on `square_topology`, of value shape
`(2,)`, with vertex values $(0, 0)$, $(1, 0)$, $(0, 1)$ and $(1, 1)$.
Using the library's `Lagrange` space, the declarations are:

```python
coordinate_space = Lagrange(square_topology, 1, (2,))
coordinates = function(coordinate_space, "coordinates")
physical_mesh = mesh(coordinates, EuclideanSpace(2))
```

- `physical_mesh.topology` is `square_topology`, with two triangle slots. `physical_mesh[0]`
  reads the first; `physical_mesh[1]` reads the second.
- The parameterisation `pull_back(coordinates, square_topology[1])` evaluates to
  $(1 - X_1, X_0 + X_1)$ for the supplied coordinates.

## Calculus on a domain

Calculus uses the Euclidean metric of the ambient coordinate domain.

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

On a coordinate domain $P$, `coordinate(P)` is the identity $X(p) = p$.

For `physical_mesh = mesh(coordinates, ambient_domain)`,
`coordinate(physical_mesh)` is defined on the physical mesh, whereas
`coordinates` is defined on its topology. Pulling the physical coordinate field
back through a cell map gives that cell's parameterisation.

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

The gradient acts componentwise. On a mesh, `Gradient`, `Divergence` and
`Curl` act cellwise, without distributional interface terms; traces and
interface integrals are explicit.

#### The frame

For a domain $P$ of dimension $k$ in $n$ ambient dimensions:

- For $k = n$, the last-axis entry `j` is the partial derivative in
  ambient coordinate $X_j$.
- For $k < n$, `Gradient(f, P)` is the tangential gradient in ambient coordinates:
  the gradient of any extension off $P$, projected orthogonally onto its
  tangent space. The result is independent of the extension.

#### Examples

Start with a scalar expression on `EuclideanSpace(3)`:

```python
coordinate_domain = EuclideanSpace(3)
coordinates = coordinate(coordinate_domain)
x, y, z = coordinates[0], coordinates[1], coordinates[2]
u = x * y * z
```

| Expression                              | Value shape | Meaning                                       |
| --------------------------------------- | ----------- | --------------------------------------------- |
| `gradient(u, coordinate_domain)`                        | `(3,)`      | The gradient $(yz, xz, xy)$, after evaluation |
| `gradient(gradient(u, coordinate_domain), coordinate_domain)`           | `(3, 3)`    | The Hessian of `u`                            |
| `gradient(coordinates, coordinate_domain)`                        | `(3, 3)`    | The identity matrix, after evaluation         |

For a scalar field `u` on a surface mesh `surface_mesh` in `EuclideanSpace(3)`,
`gradient(u, surface_mesh)` has shape `(3,)`.
`gradient(coordinate(surface_mesh), surface_mesh)` is a `(3, 3)` matrix that
orthogonally projects ambient vectors onto the surface's tangent plane.

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

### Jacobian and density

These helpers accept a `DomainMap` with a coordinate parameter domain or a
mesh `CellMap`. A topology `CellMap` has no parameterisation and raises
`TypeError`. Each helper validates its input before building the composition.
Below, `D` is the map, $P$ its $k$-dimensional parameter domain, $\psi$ its
parameterisation and $n$ its target coordinate dimension.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `jacobian(map)` | | `map: DomainMap \| CellMap` | `gradient(map.parameterisation, P)`, where `P` is `map.parameter_domain` for a `DomainMap` or `map.reference_cell` for a `CellMap` |
| `density(map)` | | `map: DomainMap \| CellMap` | `square_root(determinant(transpose(J) @ J))`, with `J = jacobian(map)` |

The Jacobian matrix $J$ has value shape `(n, k)`, and $K$ =
`inverse(jacobian(D))`, of value shape `(k, n)`, is the inverse of $J$ for
$k = n$ and the pseudo-inverse $(J^T J)^{-1} J^T$ for $k < n$. The density
$\rho$ is the factor by which $\psi$ scales the length, the area or the volume
of a small region of $P$:

$$
\rho = \sqrt{\det(J^T J)}.
$$

This equals $|\det J|$ for square $J$ and gives $1$ when $k=0$, since the
empty determinant is $1$.

$J^T J$ is `component_tensor(index_sum(J[l, a] * J[l, b], (l,)), (a, b))`, with
derived indices.

An inverse Jacobian requires $\operatorname{rank} J = k$ wherever evaluated.
The density alone does not require full rank and may vanish.

The canonical expressions for $J$, $\det J$ and $K$ are `jacobian(D)`,
`determinant(jacobian(D))` and `inverse(jacobian(D))`.

### `Integral`

`Integral(expression, integration_domain)` integrates componentwise over a
coordinate domain or mesh. The result has no spatial domain but may depend
on the integrand and the domain's geometry.

| `Integral(expression, integration_domain)` |                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------- |
| Extends          | `AbstractExpression` |
| Free function    | `integral(expression, integration_domain)`                                                                            |
| Takes            | `expression: AbstractExpression`, `integration_domain: AbstractDomain` |
| Properties | `expression`, the integrand, and `integration_domain`, the region over which it is integrated |
| Requires | `expression` has domain `integration_domain` or no domain; a mesh also accepts an expression on `integration_domain.ambient_domain` |
| Static type      | The value shape and free indices of `expression`, with no domain |
| Raises           | `ValueError` for `expression` if a requirement fails                         |

```{note}
`Integral` does not support partial integration that leaves a field on another spatial domain.
```

#### Integration over coordinate domains

On a coordinate domain $P$, `Integral(f, P)` is the Lebesgue integral
$\int_P f\,\mathrm{d}X$. In dimension zero it is the value at the single
point of $R^0$. Domain-free expressions have the
[constant-field interpretation](expressions.md#domain-free-expressions).
For `reference_triangle = UnitSimplex(2)` and `coordinates = coordinate(reference_triangle)`, `integral(coordinates[0], reference_triangle)` has
value $1/6$ and `integral(Integer(1), reference_triangle)` has value $1/2$.

#### Integration over a mesh

For a field `u` on a mesh, write

```python
integral(u, physical_mesh)
```

Users do not need to write cell pull-backs for mesh integration. Lowering
introduces them automatically; the code below explicitly writes out that expansion.

Integral lowering sums over the cell slots of `physical_mesh`. For a mesh
with one slot and a scalar field, the expansion is one reference-cell integral:

```python
cell = physical_mesh[0]
reference_field = pull_back(u, cell)
integral(reference_field * density(cell), cell.reference_cell)
```

For several cell slots, possibly of different types, the expansion is:

```python
total = Integer(0)
for cell in physical_mesh:
    reference_field = pull_back(u, cell)
    total += integral(reference_field * density(cell), cell.reference_cell)
```

The resulting `pull_back` nodes are evaluated by
[pull-back lowering](lowering.md#expanding-pull-backs).

Array-valued integrands are integrated componentwise. A domain-free array
has the same value on every cell; its indices do not select cell slots.

## Calculus through maps

For a map $\psi: P \to Q$, pull-back expresses a field on $Q$ as a field on
$P$. This lets evaluation and integration on the mapped geometry be written
using operations on the parameter domain $P$.

### Coordinate-map evaluation

Pull-back evaluates a field at the points given by a coordinate map. The
result is a field on the map's parameter domain.

For example, parameterise the diagonal edge of the unit triangle by
$\psi(t) = (1-t,t)$:

```python
reference_triangle = UnitSimplex(2)
reference_interval = UnitSimplex(1)
coordinates = coordinate(reference_triangle)
t = coordinate(reference_interval)[0]

edge_coordinates = tensor((Integer(1) - t, t))
edge = domain_map(reference_interval, edge_coordinates, reference_triangle)
u = coordinates[0] * coordinates[1]
trace = pull_back(u, edge)  # Value (1 - t) * t on reference_interval.
```

The map's target must match the field's domain: this `edge` reads fields on
`reference_triangle`. A map into `EuclideanSpace(2)` would instead read fields on that space.

### Integration through a map

To integrate a scalar field `u` through `coordinate_map = domain_map(parameter_domain, parameterisation, target_domain)`, pull it
back to the coordinate domain `parameter_domain` and multiply by the map's density:

```python
integral(pull_back(u, coordinate_map) * density(coordinate_map), parameter_domain)
```

Writing $f$, $\psi$, $P$ and $\rho$ for the field, parameterisation,
parameter domain and density, this is the integral

$$
\int_P f(\psi(X))\,\rho\,\mathrm{d}X.
$$

For the edge above, `density(edge)` is $\sqrt{2}$, so
`integral(trace * density(edge), reference_interval)` is $\sqrt{2}/6$.

Density does not depend on orientation. Circulation uses the oriented tangent
from the Jacobian; flux uses an oriented normal field (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

Integration preserves free indices in the field and the map. An explicit
`index_sum` sums the resulting integrals over a free index.

### Evaluation at points

#### Reference-cell evaluation

Let `u` be a field on `reference_cell` and `reference_position` a supplied
vector of reference coordinates. Point evaluation maps the single point of
`EuclideanSpace(0)` to that position:

```python
point_domain = EuclideanSpace(0)
point_map = domain_map(point_domain, reference_position, reference_cell)
value = integral(pull_back(u, point_map), point_domain)
```

`pull_back(u, point_map)` reads `u` at `reference_position`. Integration over
`point_domain` produces `value` with no spatial domain. The density of
`point_map` is $1$. `reference_position` is domain-free, with one component
per reference coordinate.

#### Quadrature

In finite-element assembly, quadrature is applied locally after mesh integrals
and pull-backs have been lowered. Here `reference_integrand` is scalar,
defined on `reference_cell`, and already includes the integration density.
Array-valued integrands are handled componentwise.

The supplied reference quadrature `points` have shape `(num_points, dimension)`,
where `dimension = reference_cell.dimension`; `weights` has shape
`(num_points,)`. Quadrature lowering replaces
`integral(reference_integrand, reference_cell)` with a weighted sum.
The free index `i` in `points[i]` makes `domain_map` describe a family of
point maps, one for each quadrature point:

```python
point_domain = EuclideanSpace(0)
i = index("i")
point_maps = domain_map(point_domain, points[i], reference_cell)
values = integral(pull_back(reference_integrand, point_maps), point_domain)
approximation = index_sum(weights[i] * values, (i,))
```

`approximation` equals the reference-cell integral when the quadrature rule
is exact for `reference_integrand`. Otherwise, `approximation` is a numerical
approximation to the integral.

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

A law $m$ transforms a reference value $\hat f$ into the physical value
$m_\psi(\hat f)$ using the cell map $\psi$. Let $k$ be the reference-cell
dimension and $n$ the ambient coordinate dimension:

1. Transform each trailing axis listed in `variances`: multiply by $K^T$
   for `"covariant"` or by $J$ for `"contravariant"`. Leading axes are unchanged.
2. Divide the result by $(\det J)^w$ when $k=n$, or by $\rho^w$ when
   $k<n$, where $w$ is `weight`.

Here $J$, $K$ and $\rho$ are the cell map's [Jacobian, inverse Jacobian and
density](#jacobian-and-density). Each transformed axis changes extent from
$k$ to $n$. Variance names are checked statically.

For $k=n$, the signed determinant preserves orientation information.
For $k<n$, the law uses the unsigned density, with orientation supplied separately.

Examples of transformation laws:

| Space                             | Law                                       | Mapping                |
| --------------------------------- | ----------------------------------------- | ---------------------- |
| Lagrange, vector Lagrange         | `((), 0)`                                 | identity               |
| Nédélec, $H(\mathrm{curl})$       | `(("covariant",), 0)`                     | covariant Piola        |
| Raviart–Thomas, $H(\mathrm{div})$ | `(("contravariant",), 1)`                 | contravariant Piola    |
| Discontinuous densities           | `((), 1)`                                 | $L^2$ Piola            |
| Regge                             | `(("covariant", "covariant"), 0)`         | double covariant       |
| Hellan–Herrmann–Johnson           | `(("contravariant", "contravariant"), 2)` | double contravariant   |
| $H(\mathrm{curl}\,\mathrm{div})$ matrix fields | `(("covariant", "contravariant"), 1)` | [covariant–contravariant Piola](https://defelement.org/finite-elements.html) |

The law has no orientation. A space of odd weight on a topology of lower
dimension than the ambient coordinate domain gives its own signs for each pair
of cell and facet, as the space gives the transformations of its degrees of
freedom.

#### Geometry-dependent bases

`TransformationLaw` maps field-value axes. Elements such as Hermite and Argyris
also need geometry-dependent basis transformations. The space library supplies
these through [cell evaluation](forms.md#functions-in-a-space), without requiring
a `TransformationLaw`.

### `PushForward`

`PushForward(expression, mesh, law)` carries a field from the mesh's topology
to the physical mesh. `law` specifies how its value axes and density weight
transform under each cell's coordinate map.

`PushForward` applies a `TransformationLaw` to a topology field. Geometry-dependent basis transformations are supplied through cell evaluation by the space library.

The identity law reverses pull-back cellwise; other laws also transform
values.

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
and $n$ be the topology and ambient dimensions. The pull-back of `PushForward(f, D, m)` by `cell_map(D, j)` is
$m_\psi(\hat f)$, with $\psi$ the cell and $\hat f$ =
`pull_back(f, cell_map(T, j))`. The push-forward of the coordinate function
$x$ of $D$ by `transformation_law((), 0)` is `coordinate(D)`.

The law is explicit. Libraries may supply it from a space and wrap
`push_forward`, as in the [classic FEniCS layer](../examples/forms.md#a-classic-fenics-layer-for-the-poisson-problem).
Form arguments, derivative variables and directions, and data keys remain
functions on the topology $T$.

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
| `side_domain(parent, source, *entries)` | | `parent: Mesh`, `source: AbstractTopology`, `entries: tuple[int, DomainMap]` | `mesh(pull_back(parent.coordinate_function, SideMap(parent.topology, source, entries)), parent.ambient_domain)` |

Coordinate functions may be vector-valued or tensors of scalar functions.
For a conforming interface, the library supplies parent-cell parameterisations
$x_0$, $x_1$ and reference inclusions $A_0$, $A_1$ whose physical traces agree:

$$
x_0(A_0(t)) = x_1(A_1(t))
$$

at every interval point $t$. Neither the coordinate element nor the reference
inclusions alone ensure this.

#### Traces and normals

Restriction and tangential differentiation give different gradients:

- Restricting `gradient(u, physical_mesh)` retains the full parent-cell gradient, in the
  ambient frame, as needed for a normal derivative.
- For a side map `side_relation` from a mesh `interface_mesh` to `physical_mesh`,
  `gradient(pull_back(u, side_relation), interface_mesh)` is the tangential gradient along `interface_mesh`.

#### Example

The example uses `physical_mesh` and `square_topology` of the
[`Mesh`](#mesh) example and a function `topology_field` on `square_topology`.
The diagonal has one interval slot, with an inclusion into each triangle:

```python
reference_interval = ReferenceSimplex(1)
diagonal = StaticTopology((reference_interval,), ((1, 2),))
t = coordinate(reference_interval)[0]
plus_inclusion = domain_map(reference_interval, tensor((Integer(1) - t, t)), reference_triangle)
minus_inclusion = domain_map(reference_interval, tensor((Integer(0), t)), reference_triangle)
```

Use these inclusions to read the two traces and construct the physical edge:

```python
plus = side_map(diagonal, (square_topology[0], plus_inclusion))
minus = side_map(diagonal, (square_topology[1], minus_inclusion))
plus_trace = pull_back(topology_field, plus)
minus_trace = pull_back(topology_field, minus)
jump = plus_trace - minus_trace

interface_mesh = side_domain(physical_mesh, diagonal, (0, plus_inclusion))
```

`plus` and `minus` read the two traces; `interface_mesh` takes its physical geometry from
cell 0. See [Evaluating a trace](../examples/coupling.md#evaluating-a-trace)
for their evaluation and the effect of reversing facet vertices.

### `Conormal`

`Conormal` constructs the physical unit normal to a side within its parent
cell. It maps reference conormals to vectors tangent to the parent cell and
orthogonal to the side. Changing the supplied orientation leaves the side
map and its pull-backs unchanged.

A method may instead supply a physical normal as a real vector field on the
integration domain, with one component per ambient coordinate direction.
The method supplies the orientation and normalisation. A unit normal
perpendicular to a surface in $R^3$ is supplied this way (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

| `Conormal(side_map, reference_conormals)` | |
| ------------------------------------ | - |
| Extends | `AbstractExpression` |
| Takes | `side_map: SideMap`, `reference_conormals: tuple[AbstractExpression, ...]` |
| Free function | `conormal(side_map, *reference_conormals)` |
| Requires | A mesh parent; codimension one; one real conormal per source slot, without arguments or free indices |
| Static type | Vector of the parent's ambient dimension, on `side_map.source`; no free indices or argument spaces |
| Raises | `ValueError` for an invalid count, shape, domain, scalar type or free indices |

For entry `(j, A_i)`, let `K_j` be the selected parent reference cell.
Conormal $\hat m_i$ has shape `(K_j.dimension,)` and
domain `A_i.parameter_domain` or `None`. The library supplies a nonzero vector
orthogonal to the reference inclusion with its chosen orientation. Neither
orthogonality nor nondegeneracy is established by shape validation.

Using `physical_mesh`, `interface_mesh` and `plus_inclusion` from the trace example, select the physical
parent and supply a reference conormal pointing out of its first triangle:

```python
plus_side = side_map(interface_mesh, (physical_mesh[0], plus_inclusion))
reference_conormal = tensor((Integer(1), Integer(1)))
outward = conormal(plus_side, reference_conormal)
inward = conormal(plus_side, -reference_conormal)
```

Both normals use the same evaluation map; only their orientation differs.

#### The pull-back of Conormal

Pulling `Conormal(side_relation, reference_conormals)` back to source cell `source_slot`
transforms its reference conormal by the parent cell's inverse-transpose
Jacobian, then normalises it:

```python
parent_slot, inclusion = side_relation.entries[source_slot]
parent_jacobian = jacobian(side_relation.parent[parent_slot])
inverse_at_side = pull_back(inverse(parent_jacobian), inclusion)
mapped_conormal = transpose(inverse_at_side) @ reference_conormals[source_slot]
unit_conormal = mapped_conormal / square_root(mapped_conormal @ mapped_conormal)
```

The parent Jacobian must have full column rank, and `mapped_conormal` must be nonzero.
The result lies in the parent's tangent space and is orthogonal to the mapped
facet. Its orientation follows the reference conormal.

#### The direction of the conormal

Let $k$ and $n$ be the parent mesh's topological and ambient dimensions.
Orientation follows the supplied reference conormals. For $k=n$,
normals are opposite when the two regions adjoin on opposite sides of the
interface; coincident geometry or side names alone do not imply this.
For $k<n$, matching tangent spaces are also required. On a surface, a jump
therefore uses each side's conormal, as in `u_p * n_p + u_m * n_m`.

#### Examples

Continue with the square mesh and its diagonal. For a scalar field `u` on
`physical_mesh`, the squared jump of its normal derivative, measured in the orientation
of `outward`, is:

```python
minus_side = side_map(interface_mesh, (physical_mesh[1], minus_inclusion))
field_gradient = gradient(u, physical_mesh)
grad_plus = pull_back(field_gradient, plus_side)
grad_minus = pull_back(field_gradient, minus_side)

normal_jump = (grad_plus - grad_minus) @ outward
penalty = integral(normal_jump ** Integer(2), interface_mesh)
```

Worked examples for this chapter: [an interior facet integral](../examples/coupling.md#an-interior-facet-integral), [a Bank-Weiser error estimator](../examples/estimation.md#a-bank-weiser-error-estimator), [a fluid and a solid](../examples/coupling.md#a-fluid-and-a-solid), [a 1D network in a 3D mesh](../examples/coupling.md#a-1d-network-in-a-3d-mesh), [a deformed domain](../examples/surfaces.md#a-deformed-domain).

## Differentiation with respect to functions

### `DirectionalDerivative`

`DirectionalDerivative(expression, variable, point, direction)` is the
*Gâteaux directional derivative* with a real perturbation parameter.
`point` is the evaluation value in the variable's function space;
`direction` is its perturbation.

Writing $f$, $x$, $a$ and $h$ for these four arguments, its value is

$$
\mathrm{D}(f, x, a, h) = \frac{d}{d\varepsilon} f[x := a + \varepsilon h] \Big|_{\varepsilon = 0},
$$

with $\varepsilon$ real. `variable` is bound only in `expression`.
For tuple variables, all members vary simultaneously.

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

`directional_derivative` covers differentiation with respect to parameters,
fields and geometry through the choice of variable and direction:

| Differentiation | Construction |
| --------------- | ------------ |
| Energy variation or residual linearisation | Vary the field in a supplied direction. |
| Coupled fields | Supply tuples of variables, evaluation values and directions. |
| Second and higher variations | Nest directional derivatives; see the [derivative rules](differentiation.md#directional-derivative-rules). |
| Scalar parameter derivative | Vary the parameter in the unit direction. |
| Vector or tensor derivatives | Use a family of unit directions; see [Hyperelasticity and the stress](../examples/forms.md#hyperelasticity-and-the-stress). |
| Geometry variation | Vary the coordinate function; see [Geometry preparation](differentiation.md#geometry-preparation). |
| Complex derivatives | Use real perturbations and combine directions to obtain [Wirtinger derivatives](forms.md#complex-and-wirtinger-derivatives). |

For an energy `energy` depending on a real field `u`, evaluated at `a`:

```python
variation = directional_derivative(energy, u, a, v)
linearisation = directional_derivative(variation, a, a, h)
```

Here `v` is a test direction and `h` is an increment, both independent of `a`.
Binding the directions constructs the corresponding
[derivative forms](forms.md#derivative-forms).
