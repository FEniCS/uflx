# Geometry and calculus

A mesh combines a topology with a coordinate field mapping its reference cells
into physical space. Fields can be defined on the topology or on the mesh;
pull-back and push-forward relate their reference and physical representations.
Changing the coordinate field changes the geometry while retaining the topology.

Differentiation and integration are defined on coordinate domains and meshes.
On a mesh they act cellwise, with the Jacobian relating reference and physical
derivatives and integration measures. The ambient Euclidean metric defines
tangential gradients, integration densities and unit conormals, including for
embedded curves and surfaces. Directional differentiation accounts for
dependence on fields and on the geometry itself.

## Meshes

### `Mesh`

`Mesh(coordinate_function, ambient_domain)` describes cells whose physical
coordinates are given by `coordinate_function` in `ambient_domain`. The domain
of `coordinate_function` is the topology that specifies the cell slots.

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
the mesh and preserves its requirements. Below, `D = mesh(x, Q)` and
`T = x.domain`.

#### Mesh cell maps

For a mesh `D = mesh(x, Q)` on topology `T`, `D[j]` (equivalently
`cell_map(D, j)`) selects local slot `j`. Its reference cell is `T.cell_types[j]` and its derived
parameterisation is `pull_back(x, T[j])`.

The cell map retains the selected cell when its coordinates change:
substitution reconstructs it from `D` and `j`. For a field `f` on `D`,
`pull_back(f, D[j])` is its reference-cell representation.

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

- `D.topology` is `square`, with two triangle slots. `D[0]`
  reads the first; `D[1]` reads the second.
- The parameterisation `pull_back(x, square[1])` evaluates to
  $(1 - X_1, X_0 + X_1)$ for the supplied coordinates.

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

On a coordinate domain $P$, `coordinate(P)` is the identity $X(p) = p$.

On `D = mesh(x, Q)`, `coordinate(D)` is the physical coordinate field on
`D`, whereas `x` is defined on the topology `T`. Its cell pull-back is
`pull_back(coordinate(D), D[j]) = pull_back(x, T[j])`.

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

#### Pull-backs of `Gradient`

Pull-back does not commute with `Gradient`: `PullBack(Gradient(f, Q), D)`
reads the gradient in the frame of `Q` at the mapped point.

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

For a scalar field `w` on a surface mesh `M` in `EuclideanSpace(3)`,
`gradient(w, M)` has shape `(3,)`. `gradient(coordinate(M), M)` is the
orthogonal projection onto the tangent plane, with shape `(3, 3)`.

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

$J^T J$ is `component_tensor(index_sum(J[l, a] * J[l, b], (l,)), (a, b))`, with
derived indices.

An inverse Jacobian requires $\operatorname{rank} J = k$ wherever evaluated.
The density alone does not require full rank and may vanish.

The canonical expressions for $J$, $\det J$ and $K$ are `jacobian(D)`,
`determinant(jacobian(D))` and `inverse_jacobian(D)`.

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

#### Integration over coordinate domains

On a coordinate domain $P$, `Integral(f, P)` is the Lebesgue integral
$\int_P f\,\mathrm{d}X$. In dimension zero it is the value at the single
point of $R^0$. Domain-free expressions have the
[constant-field interpretation](expressions.md#domain-free-expressions).
For `T = UnitSimplex(2)` and `X = coordinate(T)`, `integral(X[0], T)` has
value $1/6$ and `integral(Integer(1), T)` has value $1/2$.

#### Integration over mesh cells

`Integral(f, D)` is the sum of the integrals over the cell slots of `D`.
For each slot `j`, let `K_j` be its reference cell. The cell contribution is

```python
integral(pull_back(f, D[j]) * density(D[j]), K_j)
```

For the square mesh above, `integral(Integer(1), D)` is the sum of the two triangle areas, giving 1.

For a field on `D.ambient_domain`, the pull-back instead uses
`domain_map(K_j, D[j].parameterisation, D.ambient_domain)`. For such fields,
the sum equals integration over the geometric union when cell maps are
injective almost everywhere and their images overlap only on sets of zero
measure.

Array-valued integrands are integrated componentwise. A domain-free array
has the same value on every cell; its indices do not select cell slots.

## Calculus through maps

A [`DomainMap`](topology.md#domainmap) supplies a coordinate parameterisation;
[`PullBack`](topology.md#pullback) evaluates a field through it.

### Coordinate-map evaluation

Free indices in a `DomainMap` parameterisation describe a family of maps,
retaining distinct members even at identical coordinates.

To compose coordinate maps, pull back the outer parameterisation through the
inner map. Composition through a cell map keeps nested pull-backs to preserve
cell selection.

For an edge of the reference triangle `T`, choose a coordinate expression
`psi` that stretches the first coordinate by a factor of two:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
psi = tensor((Integer(2) * X[0], X[1]))
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

`i` selects samples within one cell. The cell selection determines which
trace is read at a shared point (see
[Adjacent and coincident cells](../examples/integration.md#adjacent-and-coincident-cells)).

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

A law $m$ transforms a reference value $\hat f$ into the physical value
$m_\psi(\hat f)$ using the cell map $\psi$:

1. Transform each trailing axis listed in `variances`: multiply by $K^T$
   for `"covariant"` or by $J$ for `"contravariant"`. Leading axes are unchanged.
2. Divide the result by $(\det J)^w$ for a full-dimensional cell, or by
   $\rho^w$ for an embedded cell, where $w$ is `weight`.

Here $J$, $K$ and $\rho$ are the cell map's [Jacobian, inverse Jacobian and
density](#jacobian-and-density). Each transformed axis changes extent from
the reference dimension to the ambient dimension. The determinant is signed.
Variance names are checked statically.

Examples of transformation laws:

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

The identity law reverses pull-back cellwise; other laws also transform
values. Push-forward is not a dual or adjoint operator.

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

- Restricting `gradient(f, D)` retains the full parent-cell gradient, in the
  ambient frame, as needed for a normal derivative.
- For a side map `sm` from a mesh `E` to `D`,
  `gradient(pull_back(f, sm), E)` is the tangential gradient along `E`.

#### Example

The example uses `D` and `square` of the
[`Mesh`](#mesh) example, `I = ReferenceSimplex(1)`, a
function `u` on `square`, and the diagonal of `square`:

```python
diagonal = StaticTopology((I,), ((1, 2),))
t = coordinate(I)[0]
A_plus = domain_map(I, tensor((Integer(1) - t, t)), K)
A_minus = domain_map(I, tensor((Integer(0), t)), K)
plus = side_map(diagonal, (square[0], A_plus))
minus = side_map(diagonal, (square[1], A_minus))
jump = pull_back(u, plus) - pull_back(u, minus)
E = side_domain(D, diagonal, (0, A_plus))
```

`plus` and `minus` read the two traces; `E` takes its physical geometry from
cell 0. See [Evaluating a trace](../examples/coupling.md#evaluating-a-trace)
for their evaluation and the effect of reversing facet vertices.

### Normals

A physical normal can be supplied as an ordinary real vector field on the
integration domain, with one component per ambient coordinate direction.
Fluxes and normal derivatives use this field directly; no reference conormal
is required. The supplying method determines its orientation and normalisation.
If the normal depends on geometry, that dependence must be exposed for
[geometry differentiation](differentiation.md#geometry-preparation).

### `MappedNormal`

`MappedNormal` constructs a physical unit side conormal from reference
conormals using the parent cell maps. The result is tangent to the parent cell
and orthogonal to the side. Changing the supplied orientation leaves the side
map and its pull-backs unchanged.

| `MappedNormal(side_map, reference_conormals)` | |
| ------------------------------------ | - |
| Extends | `AbstractExpression` |
| Takes | `side_map: SideMap`, `reference_conormals: tuple[AbstractExpression, ...]` |
| Free function | `mapped_normal(side_map, *reference_conormals)` |
| Requires | A mesh parent; codimension one; one real conormal per source slot, without arguments or free indices |
| Static type | Vector of the parent's ambient dimension, on `side_map.source`; no free indices or argument spaces |
| Raises | `ValueError` for an invalid count, shape, domain, scalar type or free indices |

For entry `(j, A_i)`, let `K_j` be the selected parent reference cell.
Conormal $\hat m_i$ has shape `(K_j.dimension,)` and
domain `A_i.parameter_domain` or `None`. The library supplies a nonzero vector
orthogonal to the reference inclusion with its chosen orientation. Neither
orthogonality nor nondegeneracy is established by shape validation.

Using `D`, `E` and `A_plus` from the trace example, select the physical
parent and supply a reference conormal pointing out of its first triangle:

```python
plus_mesh = side_map(E, (D[0], A_plus))
m_plus = tensor((Integer(1), Integer(1)))
outward = mapped_normal(plus_mesh, m_plus)
inward = mapped_normal(plus_mesh, -m_plus)
```

Both normals use the same evaluation map; only their orientation differs.

#### The pull-back of MappedNormal

For `sm.entries[i] = (j, A_i)`, let `B_i = cell_map(sm.parent, j)` and
$K = J^+$, the left inverse of the parent cell Jacobian. Pulling back
`MappedNormal(sm, reference_conormals)` by `cell_map(sm.source, i)` gives

$$
\frac{w}{|w|}, \qquad
w = \operatorname{pull\_back}(K^T, A_i)\,\hat m_i.
$$

This requires a full-column-rank cell Jacobian and $w\ne0$. The vector lies
in the parent's tangent space and is orthogonal to the mapped facet. The source
geometry does not choose its orientation.

#### The direction of the normal

Let $k$ and $n$ be the parent mesh's topological and ambient dimensions.
Orientation follows the supplied reference conormals. For $k=n$,
normals are opposite when the two regions adjoin on opposite sides of the
interface; coincident geometry or side names alone do not imply this.
For $k<n$, matching tangent spaces are also required. On a surface, a jump
therefore uses each side's conormal, as in `u_p * n_p + u_m * n_m`.

`MappedNormal` is a side conormal. A unit normal field to a surface in $R^3$
requires a separately supplied orientation field (see
[The flux through a surface](../examples/surfaces.md#the-flux-through-a-surface)).

#### Examples

Continue with the square mesh and its diagonal. For a scalar field `f` on
`D`, the squared jump of its normal derivative, measured in the orientation
of `outward`, is:

```python
minus_mesh = side_map(E, (D[1], A_minus))
grad_f = gradient(f, D)
i = index("i")
jump = pull_back(grad_f, plus_mesh)[i] - pull_back(grad_f, minus_mesh)[i]
penalty = integral(index_sum(jump * outward[i], (i,)) ** Integer(2), E)
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

with $\varepsilon$ real. Substitute free occurrences of $x$ in a compatible
fixed-space representation of $f$ with [geometry dependence](differentiation.md#geometry-preparation) exposed.
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

Differentiation includes dependence through geometry and spaces, with
independent field data held fixed. Moving domains require a fixed-space
representation or an explicit library derivative rule; see
[Geometry preparation](differentiation.md#geometry-preparation).

Nested derivatives are evaluated inside out. The inner derivative's point and
direction remain subject to outer differentiation; see the
[directional derivative rules](differentiation.md#directional-derivative-rules).
