# Topological meshes

UFLx describes local kernels as mathematics. A topology describes the local
cell configuration used by a kernel. Explicit side maps relate configurations.
The finite element package selects concrete mesh cells and supplies their
local coefficient and coordinate data. Mesh-wide cell counts, global cell
numbers, neighbour lookup, iteration over a mesh and global assembly are
outside the core. Runtime-sized local arrays and index maps can still describe
a complete patch calculation.

A local configuration may contain one cell, the two parents of an interface,
or the cells of a fixed macro element. Distinct local cells retain their
identities even when they have the same type or physical coordinates. A mixed
mesh uses configurations for its cell and interface types; a library can
specialise one mathematical recipe for these configurations.

Fields with geometry-independent reference-cell representations can be defined
on a topology before coordinates are supplied. A coordinate expression defines
a physical local mesh, represented by [`Mesh`](geometry.md#mesh).

```{dropdown} Spaces defined directly on a mesh
Spaces such as Hermite, Argyris and Morley have geometry-dependent
reference-cell representations and are defined directly on a physical local
mesh. The element library constructs the local basis from the coordinate map
and evaluates cell pull-backs and reference derivatives (see
[Functions in a space](forms.md#functions-in-a-space)).
```

## Topology

### `AbstractCellType`

`AbstractCellType` describes a cell type with numbered vertices.

| `AbstractCellType` | |
| ------------------ | - |
| Extends | `Base` |
| Provides | `dimension`, the nonnegative topological dimension as an `int`, and `num_vertices`, the positive number of vertices as an `int` |
| Defined by | A subclass outside the core, such as `Simplex` |

The core does not number facets or other sub-entities. A library decides which
sets of local vertices describe sub-entities and supplies their reference
maps. A vertex set identifies an entity; an ordered correspondence determines
its parameterisation. Pieces such as cuts need not be sub-entities.

A polygon's vertex count belongs to its cell type. A vertex count known only
at evaluation is not yet covered; a library specialises kernels by cell type.

### `AbstractTopology`

`AbstractTopology` describes a finite local cell configuration without physical
coordinates. Each position in `cell_types` is a distinct local cell, called a
*cell slot*. Slots identify the mathematical inputs of a kernel, not concrete
cells in the package's mesh. For reference cells, it is a labelled disjoint union; even a one-slot topology
retains its identity independently of its reusable reference cell.

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
An Alfeld macro element has three triangle slots `(K, K, K)`. A list of
alternative kernel types is managed by the library; it is not a configuration
whose cells should all be integrated together.

The number of slots is fixed when the expression is constructed. An empty
selection in the package causes no kernel calls; it does not require an empty
symbolic topology. A runtime-sized patch can instead be expressed through
packed local arrays, index tables and sums, as in the
[patch problem](../examples/patches.md#a-mean-constrained-dg-patch-correction).
The package supplies the packing and local-to-patch degree-of-freedom maps;
no runtime topology is required.

A topology whose cell types are reference cells may be the domain of an
expression or space. Its points are cell-labelled reference points; shared
vertices or sides do not identify field coefficients across slots. See
[local spaces and library bindings](forms.md#abstractfunctionspace).
A topology is not an `AbstractDomain`,
so `Coordinate`, `Gradient` and `Integral` reject it.

Identity is structural. A concrete topology stores its identifying arguments,
such as a label and the local configuration. The supplying library associates
local coefficient and coordinate inputs with each slot. The same configuration
can be evaluated repeatedly with data from different concrete cells.

## Domains

Reference cells have coordinates of their own, independently of any physical
mesh. The following domain classes provide the foundation for these reference
coordinates and for the physical domains of the next chapter.

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

A coordinate domain is a subset of $R^n$ with non-empty interior in $R^n$ whose
points are named by the coordinates of the points themselves. The core cannot
build the set from smaller objects. A library defines a new coordinate domain,
such as a reference cell, by extending this abstract class. On a coordinate
domain the coordinates are canonical: the coordinate function is the identity
(see [`Coordinate`](geometry.md#coordinate)), and the integral is the Lebesgue integral in
the coordinates (see [`Integral`](geometry.md#integral)). On a core `Mesh`,
coordinates and integration are defined through cell parameterisations.

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
It extends `AbstractCoordinateDomain` and `AbstractCellType`. The core reads no vertices of a reference cell. A library
that builds inclusions, such as from vertex positions, keeps the vertices of
its reference cells (see [The unit simplex as a reference
cell](../examples/library.md#the-unit-simplex-as-a-reference-cell)).

| `AbstractReferenceCell` |                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Extends                 | `AbstractCoordinateDomain`, `AbstractCellType` |
| Defined by              | A subclass outside the core, such as [`ReferenceSimplex`](../examples/library.md#the-unit-simplex-as-a-reference-cell)                                                   |

### `CellMap`

`CellMap(parent, cell_slot)` selects one local cell of `parent`, which is a
topology or a physical local mesh. The map goes from that cell's reference
cell into the labelled local cell of `parent`. Pulling a field back through
it gives the field's reference-cell representation on the selected slot.
The slot is part of the kernel's local configuration, not a runtime cell
number in the package's mesh.

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

A cell map of a topology has no coordinate parameterisation. For a physical
mesh its parameterisation and reconstruction are defined under
[Mesh cell maps](geometry.md#mesh-cell-maps).

| Free function | Overload | Takes | Returns |
| ------------- | -------- | ----- | ------- |
| `cell_map(parent, cell_slot)` | | `parent: AbstractTopology \| Mesh`, `cell_slot: int` | `CellMap(parent, cell_slot)` |

```python
local_cell = cell_map(parent=volume_mesh, cell_slot=0)
local_field = pull_back(field, local_cell)
# local_field has domain local_cell.reference_cell.
```

The library of the field's space supplies its local representation and data
(see [Functions in a space](forms.md#functions-in-a-space)). Ordinary indices
still describe basis functions, components, quadrature points and other local
arrays.

## Restriction

A side map is a relation from a source configuration to a parent configuration.
Each source slot selects a parent slot and a reference inclusion. The topology
itself stores no side declarations, names or normal data.

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

`E` is the shared source, with independently specified geometry. Its measure
need not be the measure induced by either parent trace. Names such as `plus`
or `"+"` are library aliases for constructed maps, not map properties.
Orientation belongs to [Normal](geometry.md#normal), not to the map identity.

### The inclusion

For `sm = SideMap(parent, source, entries)` and `entries[i] = (j, A_i)`, define
`B_i = cell_map(parent, j)`. Restriction satisfies

```python
pull_back(pull_back(f, sm), cell_map(source, i))
# becomes
pull_back(pull_back(f, B_i), A_i)
```

There is no connectivity lookup or implicit sum. The result has domain $F_i$
and the free indices of the field.

### Composition

If `inner.entries[i] = (j, A)` and `outer.entries[j] = (k, B)`, the composite
entry is `(k, C)`, where

```python
C = domain_map(A.parameter_domain, pull_back(B.parameterisation, A), B.ambient_domain)
```

The composite has `inner.source` and `outer.parent` as endpoints. Exact
intermediate endpoint identity is required; matching reference-cell types is
insufficient. Composition constructs another `SideMap`, without a new node or
an implied rule for composing outward normals. It handles subconfigurations,
permutations, repeated selections and nested inclusions uniformly.

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
domain `K`, but use different local field data. The vertex labels describe
this fixed configuration, not global mesh connectivity. The package can reuse
the kernel for many pairs of triangles.

Reference derivatives such as `gradient(u0, K)` are defined on `K`;
`gradient(u, T)` is not defined. The next chapter supplies physical local
coordinates with `D = mesh(x, Q)`.
