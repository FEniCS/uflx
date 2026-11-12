# Topological meshes

A topology describes the local cell configuration of a kernel: one cell, the
parents of an interface, or the cells of a macro element. Side maps relate
these configurations. The library selects concrete mesh cells, supplies their
coefficient and coordinate data, and assembles the results.

Reference cells have **reference geometry**, including their own coordinates.
A topology labels cell slots; when their types are reference cells, each slot
carries that reference geometry. A coordinate field maps the reference cells
into physical space, supplying **physical geometry** through
[`Mesh`](geometry.md#mesh). “Local” describes the scope of the cell configuration;
both its reference geometry and its physical geometry can be local.

Fields whose reference representations are independent of physical geometry
can be defined on a topology. Spaces whose bases depend on physical geometry
are defined directly on a mesh (see [Functions in a space](forms.md#functions-in-a-space)).

## Domains

Reference coordinates describe points within a reference cell. A physical
coordinate field maps those points into physical space. The following classes
describe the domains used for both kinds of geometry.

### `AbstractDomain`

`AbstractDomain` supplies the ambient-coordinate contract used by coordinates,
differentiation and integration. The core provides coordinate domains and `Mesh`.
Extensions may subclass `AbstractDomain`, supplying their own geometric semantics
and rewrite rules.

| `AbstractDomain` |                                                                             |
| ---------------- | --------------------------------------------------------------------------- |
| Extends          | `Base` |
| Provides         | `dimension`, the nonnegative topological dimension as an `int`, and `ambient_domain`, the `AbstractDomain` supplying coordinates for its points |
| Defined by       | `AbstractCoordinateDomain`, `Mesh`, and extension subclasses |

For an `AbstractDomain` instance `domain`, `domain.ambient_domain` supplies the coordinates
naming its points. The first coordinate domain reached by following
`ambient_domain` repeatedly is the *ambient coordinate domain*. Its dimension
is at least `domain.dimension`. A coordinate domain
is its own ambient domain.

### `AbstractCoordinateDomain`

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$,
defined by a library. Its [`Coordinate`](geometry.md#coordinate) is the identity,
and its [`Integral`](geometry.md#integral) is the Lebesgue integral in those
coordinates. On a `Mesh`, coordinates and integration are instead defined
through cell parameterisations.

| `AbstractCoordinateDomain` |                                              |
| -------------------------- | -------------------------------------------- |
| Extends                    | `AbstractDomain` |
| Provides                   | `dimension`, the nonnegative number of coordinate axes as an `int`              |
| Derived                    | `ambient_domain`, the coordinate domain itself |
| Defined by                 | A subclass, such as `EuclideanSpace` or a reference cell |

Examples:

- The unit simplex in $R^n$ is a coordinate domain. A codimension-1 face is
  not: it has empty interior in $R^n$ and is instead the image of a domain map
  from the unit simplex in $R^{n-1}$.
- A mesh of triangles is *not* a coordinate domain, because the core builds the
  mesh from a topology and a coordinate expression; its cells are parameterised
  by reference triangles.

### `EuclideanSpace`

`EuclideanSpace(dimension)` is the whole coordinate space with `dimension` axes.
Writing $n$ for `dimension`, this is $\mathbb{R}^n$.

| `EuclideanSpace(dimension)`    |                                                                 |
| ---------- | --------------------------------------------------------------- |
| Extends    | `AbstractCoordinateDomain` |
| Properties | `dimension`, the nonnegative number of coordinate axes as an `int`                                 |
| Raises     | `ValueError` if `dimension` is negative |

Worked examples for this chapter: [the unit simplex](../examples/library.md#the-unit-simplex).

## Reference cells

### `AbstractCellType`

`AbstractCellType` describes a cell type with numbered vertices.

| `AbstractCellType` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `dimension`, the nonnegative topological dimension as an `int`, and `num_vertices`, the positive number of vertices as an `int` |
| Defined by | A subclass outside the core, such as `Simplex` |

The library defines facets and other sub-entities and supplies their reference
maps. General inclusions, such as cuts, need not follow sub-entities.

### `AbstractReferenceCell`

`AbstractReferenceCell` is a cell type that is also a coordinate domain.
Libraries may supply vertex coordinates for constructing inclusions; these
are not required by the core.

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell)                                                   |

## Topology

### `AbstractTopology`

`AbstractTopology` describes a finite local cell configuration without specifying
its physical geometry. Each position in `cell_types` is a distinct local cell,
called a *cell slot*. Slots identify the mathematical inputs of a kernel.
For reference cells, the topology is a labelled disjoint union: each slot has
reference geometry and retains its identity even when reference cells or
physical coordinates coincide.

| `AbstractTopology` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `cell_types` |
| Derived | `dimension`, the dimension of the cell types |
| Defined by | A subclass outside the core |

- `cell_types` is a nonempty `tuple` of `AbstractCellType`s of the same
  dimension. Repeated types are allowed: `(K, K)` describes two distinct
  local cells of type `K`.
- For a topology `topology`, a `cell_slot` is an `int` satisfying
  `0 <= cell_slot < len(topology.cell_types)`.
  Negative numbers do not select slots from the end.

For example, `(K,)` describes a single-cell kernel, `(K, K)` the parents of
an interior-facet kernel, and `(Kt, Kq)` a triangle–quadrilateral pair.
The number of slots is fixed when the expression is constructed.

A topology of reference cells may be the domain of an expression or space.
Its points are cell-labelled reference points. The library supplies field
continuity through its [spaces and bindings](forms.md#abstractfunctionspace).
For an `AbstractTopology`, the language exposes coordinates, gradients and
integrals through its reference cells or through a mesh carrying physical geometry.

Identity is structural, determined by the concrete topology's construction
arguments, such as a label and its cell configuration. A one-slot topology is
distinct from its reference cell. The library may evaluate the same
configuration repeatedly with data from different concrete cells.

## Maps between domains

A map sends points from a starting domain to a target domain. A pull-back
reads a field on the target at those mapped points, producing a field on
the starting domain.

### `DomainMap`

`DomainMap(parameter_domain, parameterisation, ambient_domain)` describes a
geometric map by giving the target coordinates at each parameter point.
`parameter_domain` is where the map starts; `parameterisation` computes the
coordinates; `ambient_domain` is the coordinate domain containing its image.

Writing $P$, $\psi$ and $Q$ for these arguments, the map is $\psi: P \to Q$
and its geometric image is $\psi(P)$. We write $D$ for the map; its image is
always named separately. A `DomainMap` is not an
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
mesh-field evaluation uses [cell](#cellmap) or [side maps](#sidemap).

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
| [`CellMap`](#cellmap) | `reference_cell` | `parent` |
| [`SideMap`](#sidemap) | `source` | `parent` |

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

Pull-backs of terminals remain for their supplying libraries to evaluate
(see [Functions in a space](forms.md#functions-in-a-space)).

## Cell selection

### `CellMap`

`CellMap(parent, cell_slot)` selects one local cell of `parent`, which is a
topology or a mesh carrying physical geometry. The map goes from that cell's reference
cell into the labelled local cell of `parent`. Pulling a field back through
it gives the field's reference-cell representation on the selected slot.

| `CellMap(parent, cell_slot)` | |
| --------------------------- | - |
| Extends | `Base` |
| Takes | `parent: AbstractTopology \| Mesh`, `cell_slot: int` |
| Properties | `parent`, the topology or mesh carrying the field, and `cell_slot`, the local cell to read |
| Requires | `cell_slot` is a valid slot of `topology`, and `topology.cell_types[cell_slot]` is a reference cell |
| Derived | `topology`, `reference_cell`, `free_indices`, and, for a mesh, `parameterisation`, as described below |
| Raises | `ValueError` if a requirement fails |

| Derived property | Meaning |
| ---------------- | ------- |
| `topology` | `parent` itself when it is a topology; otherwise `parent.topology`. |
| `reference_cell` | `topology.cell_types[cell_slot]`, the domain of the field after pull-back. |
| `free_indices` | Always `()`; selecting a cell introduces no array index. |
| `parameterisation` | For a mesh, its coordinate expression pulled back to the selected reference cell. |

For a physical mesh, see [Mesh cell maps](geometry.md#mesh-cell-maps) for the
parameterisation and reconstruction rules.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_map(parent, cell_slot)` | | `parent: AbstractTopology \| Mesh`, `cell_slot: int` | `CellMap(parent, cell_slot)` |

Indexing a topology or mesh is shorthand for cell selection:
`parent[i]` constructs exactly `cell_map(parent, i)`. It retains the parent
and slot identity; `T.cell_types[i]` returns only the cell type. Indexing
accepts the same nonnegative integer slots as `cell_map`, without slicing
or implicit iteration.

The examples below build one configuration step by step. They use the library's
[`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell),
[`StaticTopology`](../examples/library.md#a-static-topology), and scalar
`Lagrange(topology, degree, shape)` space. All maps and pull-backs use the core
language directly.

First declare two triangle slots and a scalar field on them. Here we also
name the cell handles A and B; a library can supply these handles directly.
`function(V, "u")` declares a symbolic field in the space `V`; the library
supplies its values.

```python
K = ReferenceSimplex(2)     # Unit reference triangle.
T = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)))
V = Lagrange(T, 1, ())      # Scalar, degree-one local representation.
u = function(V, "u")        # Field on the two-slot topology T.

cell_a = T[0]               # K -> slot 0 of T.
cell_b = T[1]               # K -> slot 1 of T.
u_a = pull_back(u, cell_a)  # Cell A's field, expressed on K.
u_b = pull_back(u, cell_b)  # Cell B's field, expressed on K.
```

Both `u_a` and `u_b` have domain `K`; their cell maps select different local
field data. The vertex labels describe the correspondence along edge `(1, 2)`
in the two reference cells. The library supplies the reference inclusions for
this correspondence; physical geometry has not yet been specified.

## Restriction

A side map is a relation from a source configuration to a parent configuration.
Each source slot selects a parent slot and a reference inclusion.

### `SideMap`

| `SideMap(parent, source, entries)` | |
| -------------------------------- | - |
| Extends | `Base` |
| Takes | `parent, source: AbstractTopology \| Mesh`, `entries: tuple[tuple[int, DomainMap], ...]` |
| Properties | `parent`, `source`, `entries` |
| Derived | `free_indices`, always `()` |
| Requires | One entry per source slot; valid parent slots; inclusions without free indices and with endpoints as below |
| Raises | `ValueError` if an entry or domain condition is invalid |

For source reference cell $F_i$, entry `(j, A)` requires
`A.parameter_domain == F_i` and `A.ambient_domain == K_j`, the selected parent
reference cell. Both must be reference cells. Parent slots may repeat.
Each endpoint may independently be a topology or a mesh; for a mesh, its
topology determines its reference cells.

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `side_map(source, *entries)` | | `source: AbstractTopology \| Mesh`, `entries: tuple[CellMap, DomainMap]` | A `SideMap` with the common parent and selected slots of the cell maps |

The factory takes one `(cell, inclusion)` entry per source slot, where `cell`
is a `CellMap`. At least one entry is required, and all cell maps must have
the same parent. Missing entries or differing parents raise `ValueError`;
selectors that are not cell maps raise `TypeError`. It constructs the node above using their `cell_slot` values.
Names for the handles are ordinary Python bindings, not part of map identity.

Continue with one interval slot for the edge, ordered from vertex label 1
to vertex label 2. The library supplies `inclusion_a` and `inclusion_b`, reference inclusions
from the edge's reference interval into `K`, respecting that vertex correspondence:

```python
I = ReferenceSimplex(1)              # Unit reference interval.
S = StaticTopology((I,), ((1, 2),))  # One edge slot.

# inclusion_a, inclusion_b: library-supplied reference inclusions I -> K.
```

```{figure} ../figures/reference-slots.*
:label: reference-cell-slots
:alt: Two copies of the unit reference triangle. Cell A has vertex labels 0, 1, 2 and its diagonal edge highlighted in blue. Cell B has labels 1, 3, 2 and its left edge highlighted in amber. Both selected edges run from label 1 to label 2.
:width: 700px

These are reference triangles, not placements of cells in physical space.
Each slot has its own copy. Bold numbers are library vertex labels. The
highlighted edges show the library-supplied inclusions `inclusion_a` and `inclusion_b`;
the arrows follow the source ordering from label 1 to label 2.
```

```{mermaid}
flowchart LR
    accTitle: Two inclusions of the source interval
    accDescr: inclusion_a maps the source interval to the diagonal of the reference triangle for cell A. inclusion_b maps it to the left edge of the reference triangle for cell B.
    I["Source interval I<br/>Vertex order 1 → 2"]
    KA["K for cell A<br/>Diagonal edge"]
    KB["K for cell B<br/>Left edge"]
    I -->|inclusion_a| KA
    I -->|inclusion_b| KB
    classDef source fill:#f1f5f9,stroke:#64748b,color:#1e293b
    classDef first fill:#eff6ff,stroke:#2563a6,color:#1e293b
    classDef second fill:#fffbeb,stroke:#b45309,color:#1e293b
    class I source
    class KA first
    class KB second
    linkStyle 0 stroke:#2563a6,stroke-width:2px
    linkStyle 1 stroke:#b45309,stroke-width:2px
```

Each side map selects a parent slot and its inclusion:

```python
side_a = side_map(S, (cell_a, inclusion_a))  # Read cell A through its inclusion.
side_b = side_map(S, (cell_b, inclusion_b))  # Read cell B through its inclusion.
u_a_on_side = pull_back(u, side_a)           # Cell A's edge values, on S.
u_b_on_side = pull_back(u, side_b)           # Cell B's edge values, on S.
jump = u_a_on_side - u_b_on_side             # A field on the edge topology S.
```

The [library example](../examples/library.md#worked-edge-inclusions) shows
how to construct `inclusion_a` and `inclusion_b`.

There is one entry in each side map because `S` has one slot. In
`(cell_b, inclusion_b)`, the cell map selects the parent cell, while `inclusion_b` specifies
the reference points at which its field is read. Both inclusions respect the source vertex order,
so the two traces read corresponding points along the edge.

`side_a` and `side_b` are names for these maps. Neither chooses a normal;
normal fields and their orientation are supplied separately (see
[Normals](geometry.md#normals)).
The same construction accepts meshes carrying physical geometry as endpoints.
A source mesh then has its own physical geometry and measure, independently
of either parent trace. See [Traces on a mesh](geometry.md#traces-on-a-mesh)
for the physical-interface matching condition.

### Expressing edge values on the reference interval

To express both edge fields on the reference interval, select the sole
cell of `S`:

```python
edge = S[0]                              # I -> the edge slot of S.
u_a_on_I = pull_back(u_a_on_side, edge)  # Cell A's edge values, on I.
u_b_on_I = pull_back(u_b_on_side, edge)  # Cell B's edge values, on I.
jump_on_I = pull_back(jump, edge)        # Their difference, on I.
```

Equivalently, select the parent cell first and apply its reference inclusion:

```python
u_a_on_I = pull_back(u_a, inclusion_a)
u_b_on_I = pull_back(u_b, inclusion_b)
```

In general, for `sm.entries[i] = (j, A)`, restriction satisfies

```python
pull_back(pull_back(f, sm), cell_map(sm.source, i))
# becomes
pull_back(pull_back(f, cell_map(sm.parent, j)), A)
```

The result has domain `A.parameter_domain` and the free indices of `f`.

## Domains in the example

The domains of the resulting expressions are:

| Expression | Domain | Meaning |
| ---------- | ------ | ------- |
| `u` | `T` | Field on two labelled triangle slots |
| `u_a`, `u_b` | `K` | One selected triangle's reference representation |
| `u_a_on_side`, `u_b_on_side`, `jump` | `S` | Traces and their difference on one labelled edge slot |
| `u_a_on_I`, `u_b_on_I`, `jump_on_I` | `I` | Edge fields expressed in the interval coordinate |

Reference derivatives such as `gradient(u_a, K)` are defined on `K`.
Selecting a cell is what exposes those coordinates; `gradient(u, T)` is
not defined. A coordinate field `x` on `T` supplies physical geometry through
`D = mesh(x, Q)`, where `Q` is the physical coordinate space.
See [Geometry](geometry.md) for this construction.
