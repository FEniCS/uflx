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

`StaticTopology(cell_types, cells)` stores reference cells and library vertex
labels. Labels describe local connectivity; they do not identify field DOFs.
Side relations and orientation are constructed separately.

```python
class StaticTopology(AbstractTopology):
    identifier = "library.StaticTopology"

    def __init__(self, cell_types, cells):
        if not cell_types or len(cells) != len(cell_types):
            raise ValueError("Need one vertex tuple per local cell slot.")
        if any(K.dimension != cell_types[0].dimension for K in cell_types):
            raise ValueError("All slots must have the same dimension.")
        if any(len(v) != K.num_vertices or len(set(v)) != len(v)
               for K, v in zip(cell_types, cells)):
            raise ValueError("Need distinct vertex labels for each cell's vertices.")
        self._cell_types, self.cells = cell_types, cells

    @property
    def cell_types(self):
        return self._cell_types
```

### Reference inclusions

For affine inclusions, the library matches source vertex labels to the selected
parent cell. `affine_map` maps the first `F.dimension + 1` source vertices to
`points`; these must be affinely independent and remaining vertices must agree.

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


def topology(endpoint):
    return endpoint.topology if isinstance(endpoint, Mesh) else endpoint


def static_side_map(parent, source, *parent_slots):
    T, S = topology(parent), topology(source)
    entries = []
    for i, j in enumerate(parent_slots):
        F, K = S.cell_types[i], T.cell_types[j]
        positions = tensor(tuple(K.vertices[T.cells[j].index(v)] for v in S.cells[i]))
        entries.append((parent[j], domain_map(F, affine_map(F, positions), K)))
    return side_map(source, *entries)
```

Examples may collect reference maps in a dictionary. This convenience is outside
the core; names do not enter the identity of a topology or map:

```python
def static_sides(source, selections):
    # selections[i] lists (name, parent_topology, parent_slot) for source slot i.
    return {
        name: static_side_map(parent, source, *(
            next(j for n, T, j in row if n == name and T == parent)
            for row in selections
        ))
        for name, parent, _ in selections[0]
    }
```

#### Worked edge inclusions

The [topological-mesh walkthrough](../language/topology.md#sidemap) takes
`inclusion_a` and `inclusion_b` as library inputs. For its two unit reference triangles, the
library can construct them explicitly:

```python
K = ReferenceSimplex(2)
I = ReferenceSimplex(1)
t = coordinate(I)[0]

# Source vertex order is 1 -> 2 in each parent cell.
inclusion_a = domain_map(I, tensor((Integer(1) - t, t)), K)  # Cell A: diagonal edge.
inclusion_b = domain_map(I, tensor((Integer(0), t)), K)      # Cell B: left edge.
```

At `t = 0`, `inclusion_a` reads `(1, 0)` and `inclusion_b` reads `(0, 0)` in their respective
parent cells; both points have vertex label 1. At `t = 1`, both read `(0, 1)`,
with label 2. The library's `static_side_map` constructs these inclusions from
the vertex correspondence rather than requiring the caller to write formulae.

For local fields $u_a(\xi,\eta)=\xi+2\eta$ and
$u_b(\xi,\eta)=3\xi-\eta$, these inclusions give traces $1+t$ and $-t$,
and a jump $1+2t$.

### Reference orientation

For a planar facet of a convex reference cell, this library chooses the
conormal away from the parent centroid. Interior cuts require a different
orientation policy. The core only transports and normalises the supplied vector.

```python
def outward_conormals(relation):
    result = []
    for _, A in relation.entries:
        F, K = A.parameter_domain, A.ambient_domain
        J = jacobian(A)
        m = tensor(tuple(
            determinant(tensor(tuple(
                tensor((Integer(1 if q == r else 0),
                        *(J[q, c] for c in range(F.dimension))))
                for q in range(K.dimension)
            )))
            for r in range(K.dimension)
        ))
        centre = K.vertices[0]
        for v in range(1, K.num_vertices):
            centre = centre + K.vertices[v]
        centre = centre / Integer(K.num_vertices)
        away = m @ (A.parameterisation - centre)
        result.append((away / absolute(away)) * m)
    return tuple(result)


def outward_normal(relation):
    return conormal(relation, *outward_conormals(relation))
```

These routines assume nondegenerate affine facet inclusions. They are library
policies, not requirements on general core relations.

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
reference-inclusion constructors read it (see [Reference inclusions](#reference-inclusions)).

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


def restrict(obj, E, relation):
    # The legacy apply_restrictions, at construction: the largest subtrees on the
    # mesh are pulled back to E, and the facet normal of the mesh becomes the
    # normal of E.
    if is_facet_normal(obj):
        return outward_normal(relation)
    if isinstance(obj, tuple):
        return tuple(restrict(item, E, relation) for item in obj)
    if not isinstance(obj, AbstractExpression) or obj.domain is None:
        return obj
    if not contains(obj, is_facet_normal):
        return pull_back(obj, relation)
    return obj.__class__(*restrict(obj.construction_arguments, E, relation))


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
        reference = None
        if self.subdomain_id != "everywhere":
            reference = self.subdomain_data[self.subdomain_id]
        elif self.integral_type == "ds":
            reference = mesh.boundary
        if reference is not None:
            E = side_domain(mesh.domain, reference.source, *reference.entries)
            relation = SideMap(mesh.domain, E, reference.entries)
            integrand = restrict(integrand, E, relation)
        else:
            E = integrand.domain if mesh is None else mesh.domain
        expression = integral(integrand, E)
        args = arguments(expression)
        if not args:
            return expression
        # Legacy forms are antilinear in the test function.
        return form(expression, (conjugate(args[0]),) + args[1:])


dx = Measure("dx")
ds = Measure("ds")
```
