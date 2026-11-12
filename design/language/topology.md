# Topological meshes

A topology describes the local cell configuration of a kernel: one cell, the
parents of an interface, or the cells of a macro element. Side maps relate
these configurations. The library selects concrete mesh cells, supplies their
coefficient and coordinate data, and assembles the results.

Fields with geometry-independent reference representations can be defined on
a topology. A coordinate expression gives it a physical geometry, represented
by [`Mesh`](geometry.md#mesh). Spaces with geometry-dependent bases are defined
directly on a mesh (see [Functions in a space](forms.md#functions-in-a-space)).

## Topology

### `AbstractCellType`

`AbstractCellType` describes a cell type with numbered vertices.

| `AbstractCellType` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `dimension`, the nonnegative topological dimension as an `int`, and `num_vertices`, the positive number of vertices as an `int` |
| Defined by | A subclass outside the core, such as `Simplex` |

The library defines facets and other sub-entities and supplies their reference
maps. General inclusions, such as cuts, need not follow sub-entities.

### `AbstractTopology`

`AbstractTopology` describes a finite local cell configuration without physical
coordinates. Each position in `cell_types` is a distinct local cell, called a
*cell slot*. Slots identify the mathematical inputs of a kernel. For reference
cells, the topology is a labelled disjoint union: each slot retains its identity even
when reference cells or physical coordinates coincide.

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
A topology is not an `AbstractDomain`: coordinates, gradients and integrals
are defined on its reference cells or on a physical mesh.

Identity is structural, determined by the concrete topology's construction
arguments, such as a label and its cell configuration. A one-slot topology is
distinct from its reference cell. The library may evaluate the same
configuration repeatedly with data from different concrete cells.

## Domains

Reference cells have coordinates independently of any physical mesh. The
following classes describe reference and physical domains.

### `AbstractDomain`

`AbstractDomain` represents a domain, the set on which an expression is defined.
The core provides coordinate domains and `Mesh`. Extensions may subclass
`AbstractDomain`, supplying their own geometric semantics and rewrite rules.

| `AbstractDomain` |                                                                             |
| ---------------- | --------------------------------------------------------------------------- |
| Extends          | `Base` |
| Provides         | `dimension`, the nonnegative topological dimension as an `int`, and `ambient_domain`, the `AbstractDomain` supplying coordinates for its points |
| Defined by       | `AbstractCoordinateDomain`, `Mesh`, and extension subclasses |

For a domain `domain`, `domain.ambient_domain` supplies the coordinates
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

### `AbstractReferenceCell`

`AbstractReferenceCell` is a cell type that is also a coordinate domain.
Libraries may supply vertex coordinates for constructing inclusions; these
are not required by the core.

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell)                                                   |

### `CellMap`

`CellMap(parent, cell_slot)` selects one local cell of `parent`, which is a
topology or a physical local mesh. The map goes from that cell's reference
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

```python
local_cell = cell_map(parent=volume_mesh, cell_slot=0)
local_field = pull_back(field, local_cell)
# local_field has domain local_cell.reference_cell.
```

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
| `side_map(parent, source, *entries)` | | `parent, source: AbstractTopology \| Mesh`, `entries: tuple[int, DomainMap]` | `SideMap(parent, source, entries)` |
| `compose(outer, inner)` | | `outer, inner: SideMap` | The composed `SideMap`, requiring `inner.parent == outer.source` |

The variadic factory collects entries into the node's tuple:

```python
plus = side_map(D, E, (0, A_plus))
minus = side_map(D, E, (1, A_minus))
jump = pull_back(u, plus) - pull_back(u, minus)
```

`E` has its own geometry and measure, independently of either parent trace.
`plus` and `minus` name the constructed maps. Orientation is supplied
separately to [Normal](geometry.md#normal).

### The inclusion

For `sm = SideMap(parent, source, entries)` and `entries[i] = (j, A_i)`, define
`B_i = cell_map(parent, j)`. Restriction satisfies

```python
pull_back(pull_back(f, sm), cell_map(source, i))
# becomes
pull_back(pull_back(f, B_i), A_i)
```

The result has domain $F_i$ and the free indices of the field.

### Composition

If `inner.entries[i] = (j, A)` and `outer.entries[j] = (k, B)`, the composite
entry is `(k, C)`, where

```python
C = domain_map(A.parameter_domain, pull_back(B.parameterisation, A), B.ambient_domain)
```

The composite is a `SideMap` from `inner.source` to `outer.parent`, requiring
`inner.parent == outer.source`. This covers permutations, repeated selections
and nested inclusions.

## A field on a topological mesh

Use the library's [static topology](../examples/library.md#a-static-topology)
to describe two distinct triangle slots:

```python
K = ReferenceSimplex(2)
T = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)))
u0 = pull_back(u, cell_map(T, 0))
u1 = pull_back(u, cell_map(T, 1))
```

Here `u` is a scalar function in a library space on `T`. Both expressions have
domain `K`, but use different local field data. The library can reuse the
kernel for many pairs of triangles.

Reference derivatives such as `gradient(u0, K)` are defined on `K`;
`gradient(u, T)` is not defined. The next chapter supplies physical local
coordinates with `D = mesh(x, Q)`.
