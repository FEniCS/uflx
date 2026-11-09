# Reference library for the examples

These supporting classes are outside the core. The simplex cell type,
static topology and reference-simplex definitions carry the original draft's
LLM-generated-content designation; the centroid example was also LLM-proposed.

## Transformation laws

The examples use the library overload `transformation_law(V)` to obtain a
law for a supported space. It derives the law from the space's definition
and raises `TypeError` for an unsupported space. The two-argument form
`transformation_law(variances, weight)` is the core constructor.

## The unit simplex

Define the unit simplex in $R^n$,

$\{p \in R^n : p_i \ge 0 \text{ for every } i, \ \sum_i p_i \le 1\}$,

as a concrete `AbstractCoordinateDomain`.

```python
class UnitSimplex(AbstractCoordinateDomain):
    identifier = "library.UnitSimplex"

    def __init__(self, n: int):
        if type(n) is not int:
            raise TypeError("Dimension must be an int.")
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n
```

`UnitSimplex(1)` is the unit interval $[0, 1]$. `UnitSimplex(2)` is the unit
triangle with vertices $(0, 0)$, $(1, 0)$ and $(0, 1)$.

## The simplex as a cell type

`Simplex(n)` is the simplex of dimension $n \ge 0$ as a cell type, with
$n + 1$ vertices. `Simplex(0)` is a point.

```python
class Simplex(AbstractCellType):
    identifier = "library.Simplex"

    def __init__(self, n: int):
        if type(n) is not int:
            raise TypeError("Dimension must be an int.")
        if n < 0:
            raise ValueError("Dimension must be non-negative.")
        self._n = n

    @property
    def dimension(self) -> int:
        return self._n

    @property
    def num_vertices(self) -> int:
        return self._n + 1
```

## A static topology

`StaticTopology(cell_types, cells, sides)` describes a fixed local configuration,
mainly for examples and macro elements. `cells[j][v]` is the label of local
vertex `v` of slot `j`. These labels identify vertices within the example;
they are not mesh-wide cell numbers or runtime connectivity.

A side `(name, T, j)` explicitly names the parent slot. The library pairs the
side cell's vertex labels with that parent's vertex labels and builds an
affine map through the corresponding reference vertices. This helper supports
whole cells and sub-entities with nondegenerate affine inclusions; general
cuts use another library implementation of `side_data`. Reference normals
below assume planar facets of convex reference cells.

```python
class StaticTopology(AbstractTopology):
    identifier = "library.StaticTopology"

    def __init__(
        self,
        cell_types: tuple[AbstractCellType, ...],
        cells: tuple[tuple[int, ...], ...],
        sides: tuple[tuple[tuple[str, AbstractTopology, int], ...], ...],
    ):
        if not cell_types or len(cells) != len(cell_types):
            raise ValueError("Need one vertex tuple per local cell slot.")
        if any(K.dimension != cell_types[0].dimension for K in cell_types):
            raise ValueError("All slots must have the same dimension.")
        if any(
            len(v) != K.num_vertices or len(set(v)) != len(v)
            for K, v in zip(cell_types, cells)
        ):
            raise ValueError("Need distinct vertex labels for each cell's vertices.")
        if sides and len(sides) != len(cell_types):
            raise ValueError("Need one side tuple per local cell slot.")
        for i, entries in enumerate(sides):
            if len({name for name, _, _ in entries}) != len(entries):
                raise ValueError("Side names must be distinct within a slot.")
            for _, T, j in entries:
                if not 0 <= j < len(T.cell_types):
                    raise ValueError("Invalid parent slot.")
                if not isinstance(
                    cell_types[i], AbstractReferenceCell
                ) or not isinstance(T.cell_types[j], AbstractReferenceCell):
                    raise ValueError("Sides require reference cells.")
                if not set(cells[i]) <= set(T._cells[j]):
                    raise ValueError(
                        "Side vertices must belong to the declared parent."
                    )
        self._cell_types = cell_types
        self._cells = cells
        self._sides = sides

    @property
    def cell_types(self) -> tuple[AbstractCellType, ...]:
        return self._cell_types

    @property
    def sides(self) -> tuple[tuple[tuple[str, AbstractTopology, int], ...], ...]:
        return self._sides

    def side_data(self, cell_slot: int, side: str):
        return static_side_data(self, cell_slot, side)
```

With `K = ReferenceSimplex(2)` and `I = ReferenceSimplex(1)`
(defined below):

```python
# Two triangles forming a square.
square = StaticTopology((K, K), ((0, 1, 2), (1, 3, 2)), ())

# Three triangles meeting at vertex 3 (an Alfeld split).
alfeld = StaticTopology((K, K, K), ((0, 1, 3), (1, 2, 3), (2, 0, 3)), ())

# Four boundary facets, each seen from one triangle.
boundary = StaticTopology(
    (I, I, I, I),
    ((0, 1), (1, 3), (3, 2), (2, 0)),
    tuple((("inside", square, j),) for j in (0, 1, 1, 0)),
)

# One diagonal facet, seen from two triangles.
diagonal = StaticTopology((I,), ((1, 2),), ((("+", square, 0), ("-", square, 1)),))
```

The four configurations have two, three, four and one local slots,
respectively. The parent slot is explicit in each side declaration; no
neighbour search or table read occurs inside the symbolic kernel.

### Side data of a static topology

The library builds [side data](../language/topology.md#side-data) from reference
vertices, which the core does not read. For a $k$-dimensional reference cell
$F$ with vertices $V_v$, `affine_map(F, points)` maps the first $k+1$ vertices
to rows $p_v$ of `points`:

$$
X \mapsto p_0 + \sum_{d < k} (p_{d+1} - p_0) \, \big(E^{-1} (X - V_0)\big)_d,
$$

Here $E$ has columns $V_{d+1}-V_0$. The implementation uses a derived index:

```python
def affine_map(F, points):
    l = Index(("library.affine_map", F, points, 0))
    k, V = F.dimension, F.vertices
    if k == 0:
        return component_tensor(points[0, l], (l,))
    E = tensor(tuple(tensor(tuple(V[d + 1, r] - V[0, r] for d in range(k))) for r in range(k)))
    Y = inverse(E) @ (coordinate(F) - V[0])
    total = points[0, l]
    for d in range(k):
        total = total + (points[d + 1, l] - points[0, l]) * Y[d]
    return component_tensor(total, (l,))
```

Those $k+1$ vertices must be affinely independent; later rows of `points`
are ignored. A sub-entity with vertex numbers $(w_0,w_1,\dots)$ in `K` uses
row `K.vertices[w_v]` for point $v$. For example:

| Parent | Sub-entity vertices | Reference cell | Inclusion |
| ------ | ------------------- | -------------- | --------- |
| Unit triangle | `(1, 2)` | Unit interval | $t\mapsto(1-t,t)$ |
| Unit interval | `(1,)` | Reference point | The endpoint $1$ |

The following helper constructs the pair `(inclusion, reference_normal)`.
It uses the same vertex correspondence on both sides, so both traces are
parameterised by the same reference point. The reference normal need not
have unit length.

```python
def static_side_data(S, i, side):
    if not 0 <= i < len(S.cell_types):
        raise ValueError("Invalid local slot.")
    entries = S.sides[i] if S.sides else ()
    matches = [(T, j) for name, T, j in entries if name == side]
    if len(matches) != 1:
        raise ValueError("Unknown side.")
    T, j = matches[0]
    F, K = S.cell_types[i], T.cell_types[j]
    parent_vertices = T._cells[j]
    points = tensor(tuple(K.vertices[parent_vertices.index(v)] for v in S._cells[i]))
    inclusion = domain_map(F, affine_map(F, points), K)
    if F.dimension != K.dimension - 1:
        return inclusion, None

    k = K.dimension
    J = jacobian(inclusion)
    cofactors = []
    for r in range(k):
        C = tensor(
            tuple(
                tensor(
                    (Integer(1 if q == r else 0),)
                    + tuple(J[q, d] for d in range(k - 1))
                )
                for q in range(k)
            )
        )
        cofactors.append(determinant(C))
    m0 = tensor(tuple(cofactors))
    centroid = K.vertices[0]
    for v in range(1, K.num_vertices):
        centroid = centroid + K.vertices[v]
    centroid = centroid / Integer(K.num_vertices)
    away = m0 @ (points[0] - centroid)
    return inclusion, away / absolute(away) * m0
```

This formula assumes full-rank affine maps onto boundary facets of convex
parents, so `away` is nonzero. Interior cuts need a supplier-chosen orientation.
For nonsimplex cells the first `F.dimension + 1` vertices must be affinely
independent and the remaining vertices must agree with the resulting map.
These are assumptions of this example helper, not limitations on core side data.

## The unit simplex as a reference cell

`ReferenceSimplex(n)` is the cell type `Simplex(n)` as a reference cell, the
unit simplex in $R^n$. Vertex 0 is the origin, and vertex $v \ge 1$ is the unit
vector $e_{v-1}$. `ReferenceSimplex(0)` is the reference point.

```python
def simplex_vertices(n: int) -> AbstractExpression:
    if type(n) is not int:
        raise TypeError("Dimension must be an int.")
    if n < 0:
        raise ValueError("Dimension must be non-negative.")
    if n == 0:
        return Zero((1, 0), (), (), real)
    return tensor(
        tuple(
            tensor(tuple(Integer(1 if k == v - 1 else 0) for k in range(n))) for v in range(n + 1)
        )
    )


class ReferenceSimplex(Simplex, AbstractReferenceCell):
    identifier = "library.ReferenceSimplex"

    def __init__(self, n: int):
        super().__init__(n)
        self._vertices = simplex_vertices(n)

    @property
    def vertices(self) -> AbstractExpression:
        return self._vertices
```

`vertices` is a property of the library, not of `AbstractReferenceCell`. The
side data of a `StaticTopology` read it (see [Side data of a static
topology](#side-data-of-a-static-topology)).

`ReferenceSimplex(2)` and `Simplex(2)` are different cell types, so a topology
of `Simplex(2)` cells has no reference cells. `ReferenceSimplex(2)` and
`UnitSimplex(2)` are the same set but different coordinate domains, so an
expression on one is not an expression on the other. Likewise,
`ReferenceSimplex(0)` and `EuclideanSpace(0)` are different coordinate domains.

## Implementation of the FEniCS layer

The [Poisson example](forms.md#a-classic-fenics-layer-for-the-poisson-problem) uses
these wrappers. `Lagrange(topology, degree, shape)`, `Constants(shape)` and
`Fields(domain, shape)` stand for library spaces: finite element functions,
domainless constants and fields without an element, respectively. The plain
Python wrapper classes are outside the core.

```python
# uflx/legacy.py
from uflx import (
    AbstractExpression,
    Base,
    Function,
    conjugate,
    curl as core_curl,
    divergence,
    form,
    function,
    gradient,
    index,
    index_sum,
    integral,
    normal,
    pull_back,
    push_forward,
    side_domain,
    side_map,
)


class FunctionSpace:
    # A space on the topology of the mesh, with the mesh that its functions are pushed to.
    def __init__(self, mesh, family, degree, shape=()):
        self.mesh = mesh
        self.space = Lagrange(mesh.domain.topology, degree, shape)


def VectorFunctionSpace(mesh, family, degree):
    return FunctionSpace(mesh, family, degree, (mesh.domain.ambient_domain.dimension,))


def TestFunction(V):
    return push_forward(
        function(V.space, ("argument", 0)), V.mesh.domain, transformation_law(V.space)
    )


def TrialFunction(V):
    return push_forward(
        function(V.space, ("argument", 1)), V.mesh.domain, transformation_law(V.space)
    )


def Coefficient(V, name):
    return push_forward(function(V.space, name), V.mesh.domain, transformation_law(V.space))


def Constant(name, shape=()):
    return function(Constants(shape), name)


def FacetNormal(mesh):
    # A placeholder on the mesh, which restrict replaces by the normal of the facets.
    n = mesh.domain.ambient_domain.dimension
    return function(Fields(mesh.domain, (n,)), "facet normal")


def is_facet_normal(obj):
    return isinstance(obj, Function) and obj.label == "facet normal"


def grad(f):
    return gradient(f, f.domain)


def div(w):
    return divergence(w, w.domain)


def curl(w):
    return core_curl(w, w.domain)


def inner(a, b):
    # For scalars and vectors, conjugating the second operand as legacy UFL does.
    if a.value_shape == ():
        return a * conjugate(b)
    i = index("inner")
    return index_sum(a[i] * conjugate(b[i]), (i,))


def contains(obj, predicate):
    if predicate(obj):
        return True
    if isinstance(obj, tuple):
        return any(contains(item, predicate) for item in obj)
    return isinstance(obj, Base) and any(
        contains(item, predicate) for item in obj.construction_arguments
    )


def arguments(expression):
    # The legacy extract_arguments: the functions labelled ("argument", n), by n,
    # which the walk finds inside the push-forwards.
    found = {}

    def visit(obj):
        if isinstance(obj, Function) and isinstance(obj.label, tuple):
            if obj.label[0] == "argument":
                found[obj.label[1]] = obj
        if isinstance(obj, tuple):
            for item in obj:
                visit(item)
        elif isinstance(obj, Base):
            for item in obj.construction_arguments:
                visit(item)

    visit(expression)
    return tuple(found[n] for n in sorted(found))


def restrict(obj, E):
    # The legacy apply_restrictions, at construction: the largest subtrees on the
    # mesh are pulled back to E, and the facet normal of the mesh becomes the
    # normal of E.
    if is_facet_normal(obj):
        return normal(side_map(obj.domain, E, "inside"))
    if isinstance(obj, tuple):
        return tuple(restrict(item, E) for item in obj)
    if not isinstance(obj, AbstractExpression) or obj.domain is None:
        return obj
    if not contains(obj, is_facet_normal):
        return pull_back(obj, side_map(obj.domain, E, "inside"))
    return obj.__class__(*restrict(obj.construction_arguments, E))


class Mesh:
    def __init__(self, domain, boundary):
        self.domain, self.boundary = domain, boundary


class Measure:
    def __init__(self, integral_type, domain=None, subdomain_id="everywhere", subdomain_data=None):
        self.integral_type = integral_type
        self.domain = domain
        self.subdomain_id = subdomain_id
        self.subdomain_data = subdomain_data

    def __call__(self, subdomain_id="everywhere", domain=None, subdomain_data=None):
        # dx(1), or dx(domain=mesh, subdomain_data=markers), fills the slots.
        return Measure(
            self.integral_type,
            self.domain if domain is None else domain,
            subdomain_id,
            self.subdomain_data if subdomain_data is None else subdomain_data,
        )

    def __rmul__(self, integrand):
        # expr * dx is a form in legacy UFL, so the measure binds the arguments.
        mesh = self.domain
        if self.subdomain_id != "everywhere":
            # The marker value names a topology with a side "inside" into the mesh.
            E = side_domain(mesh.domain, self.subdomain_data[self.subdomain_id], "inside")
        elif self.integral_type == "ds":
            E = side_domain(mesh.domain, mesh.boundary, "inside")
        else:
            E = integrand.domain if mesh is None else mesh.domain
        if E != integrand.domain:
            integrand = restrict(integrand, E)
        expression = integral(integrand, E)
        args = arguments(expression)
        if not args:
            return expression
        # Legacy forms are antilinear in the test function.
        return form(expression, (conjugate(args[0]),) + args[1:])


dx = Measure("dx")
ds = Measure("ds")
```
