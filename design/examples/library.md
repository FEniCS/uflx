# Reference library for the examples

These supporting classes are outside the core. The simplex cell type,
static topology and reference-simplex definitions carry the original draft's
LLM-generated-content designation; the centroid example was also LLM-proposed.

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

`StaticTopology(cell_types, cells, sides)` stores fixed connectivity, mainly
for testing UFLx ideas, including Alfeld splits. It has no coordinates and
does not check its sides; vertex positions come from a mesh's coordinate
function.

| Convention | Meaning |
| ---------- | ------- |
| `cells[j][c][i]` | Topology vertex number of local vertex `i` in cell `c` of type `j`. |
| Vertex pairing | Each side uses the same vertex numbers as its target topology. |
| Inclusions | Affine maps through the paired reference vertices; sides cover sub-entities, not pieces of cells. |
| Two cells in the same group | The lower-numbered cell is on the first side. |
| Two cells in different groups | Each side's group selects its cell. |
| Sub-entities below facet dimension | Select the lowest-numbered containing cell. |

The library derives the containing cells, inclusions and facet normals as
described in [Side tables](#side-tables-of-a-static-topology).

```python
class StaticTopology(AbstractTopology):
    identifier = "library.StaticTopology"

    def __init__(
        self,
        cell_types: tuple[AbstractCellType, ...],
        cells: tuple[tuple[tuple[int, ...], ...], ...],
        sides: tuple[tuple[tuple[str, AbstractTopology, int], ...], ...],
    ):
        if len(cells) != len(cell_types):
            raise ValueError("Need one tuple of cells for each cell type.")
        self._cell_types = cell_types
        self._cells = cells
        self._sides = sides

    @property
    def cell_types(self) -> tuple[AbstractCellType, ...]:
        return self._cell_types

    @property
    def num_cells(self) -> tuple[int, ...]:
        return tuple(len(c) for c in self._cells)

    @property
    def sides(self) -> tuple[tuple[tuple[str, AbstractTopology, int], ...], ...]:
        return self._sides
```

With `K = Simplex(2)` and `I = Simplex(1)`:

```python
# Two triangles forming a square.
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())

# Three triangles meeting at vertex 3 (an Alfeld split).
alfeld = StaticTopology((K,), (((0, 1, 3), (1, 2, 3), (2, 0, 3)),), ())

# Four boundary facets, each seen from one triangle.
boundary = StaticTopology(
    (I,), (((0, 1), (1, 3), (3, 2), (2, 0)),), ((("inside", square, 0),),)
)

# One diagonal facet, seen from two triangles.
diagonal = StaticTopology(
    (I,), (((1, 2),),), ((("+", square, 0), ("-", square, 0)),)
)
```

Facets are matched by vertex sets, using the facets of `K`. The square's
diagonal and the Alfeld split's three facets containing vertex 3 each belong
to two cells; their remaining facets belong to one. The four topologies have
`num_cells` `(2,)`, `(3,)`, `(4,)` and `(1,)`, respectively.

### Side tables of a static topology

The library builds [side tables](../language.md#side-tables) from reference
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

For group $i$ with reference cell `F_i`, let `e` be its cell index and
`positions[e, v]` the position of vertex `v` in the side's reference cell
`K_j`. With derived indices `v` and `k`, its inclusion is:

```python
points = component_tensor(positions[e, v, k], (v, k))
inclusion = domain_map(F_i, affine_map(F_i, points), K_j)
```

For a group of facets, the reference normal is the cofactor vector of the
Jacobian of the inclusion, turned away from the centroid of $K_j$:

```python
J_A = jacobian(inclusion)  # value shape (k, k - 1)
C_r = tensor(tuple(
    tensor((Integer(1 if q == r else 0), J_A[q, 0], ..., J_A[q, k - 2]))
    for q in range(k)
))
m0 = tensor((determinant(C_0), ..., determinant(C_(k - 1))))
centroid = (V[0] + V[1] + ... + V[m_j - 1]) / Integer(m_j)
away = m0 @ (points[0] - centroid)
reference_normal = away / absolute(away) * m0
```

Here `V = K_j.vertices`, $m_j$ is its vertex count, and the sums nest to the
left. All sides preserve the vertex order of cell `e`, so their inclusions
parameterise it in the same way.

## The unit simplex as a reference cell

`ReferenceSimplex(n)` is the cell type `Simplex(n)` as a reference cell, the
unit simplex in $R^n$. Vertex 0 is the origin, and vertex $v \ge 1$ is the unit
vector $e_{v-1}$. `ReferenceSimplex(0)` is the reference point.

```python
class ReferenceSimplex(Simplex, AbstractReferenceCell):
    identifier = "library.ReferenceSimplex"

    @property
    def vertices(self) -> AbstractExpression:
        n = self.dimension
        if n == 0:
            return Zero((1, 0), (), ())
        return tensor(
            tuple(
                tensor(tuple(Integer(1 if k == v - 1 else 0) for k in range(n)))
                for v in range(n + 1)
            )
        )
```

`vertices` is a property of the library, not of `AbstractReferenceCell`. The
side tables of a `StaticTopology` read it (see [Side tables of a static
topology](#side-tables-of-a-static-topology)).

`ReferenceSimplex(2)` and `Simplex(2)` are different cell types, so a topology
of `Simplex(2)` cells has no reference cells. `ReferenceSimplex(2)` and
`UnitSimplex(2)` are the same set but different coordinate domains, so an
expression on one is not an expression on the other. Likewise,
`ReferenceSimplex(0)` and `RD(0)` are different coordinate domains.

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
    return push_forward(function(V.space, ("argument", 0)), V.mesh.domain)


def TrialFunction(V):
    return push_forward(function(V.space, ("argument", 1)), V.mesh.domain)


def Coefficient(V, name):
    return push_forward(function(V.space, name), V.mesh.domain)


def Constant(name, shape=()):
    return function(Constants(shape), name)


def FacetNormal(mesh):
    # A placeholder on the mesh, which restrict replaces by the normal of the facets.
    n = mesh.domain.ambient_domain.dimension
    return function(Fields(mesh.domain, (n,)), "facet normal")


def is_facet_normal(obj):
    return isinstance(obj, Function) and obj.label == "facet normal"


def grad(f):
    return gradient(f, f.support)


def div(w):
    return divergence(w, w.support)


def curl(w):
    return core_curl(w, w.support)


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
        return normal(side_map(obj.support, E, "inside"))
    if isinstance(obj, tuple):
        return tuple(restrict(item, E) for item in obj)
    if not isinstance(obj, AbstractExpression) or obj.support is None:
        return obj
    if not contains(obj, is_facet_normal):
        return pull_back(obj, side_map(obj.support, E, "inside"))
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
            E = side_domain(mesh.domain, self.subdomain_data[self.subdomain_id])
        elif self.integral_type == "ds":
            E = side_domain(mesh.domain, mesh.boundary)
        else:
            E = integrand.support if mesh is None else mesh.domain
        if E != integrand.support:
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
