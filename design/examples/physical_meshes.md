# Physical polynomial spaces on polygon meshes

A library can define $V_h=\{u:\ u|_K\in\mathbb P_p(K)\}$ directly on physical
polygons, with independent polynomial coefficients on each cell. No reference
polygon is needed. A subdivision used for quadrature does not define the space.

Here a triangle and a quadrilateral share the edge between vertices 1 and 2.
Each cell lists its vertex numbers in counterclockwise order. The illustrative
`polygon_topology` factory describes two local slots, with a triangle and
a quadrilateral cell type. Mesh, space and
physical calculus functions below are also extension functions; algebra and
form binding use core names:

```python
vertices = ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1))
topology = polygon_topology(  # extends AbstractTopology
    cells=((0, 1, 2), (1, 3, 4, 2)),
)
D = polygon_mesh(topology, vertices)  # extends AbstractDomain
V = physical_polynomial_space(D, degree=1)  # extends AbstractFunctionSpace
u, v = function(V, "u"), function(V, "v")

i = index("i")
a = form(
    integral(index_sum(gradient(v, D)[i] * gradient(u, D)[i], (i,)), D),
    (v, u),
)
```

A language extension can describe physical integration, derivatives and
traces, reusing core algebra and form binding where their contracts apply.
Its compiler can lower volume terms directly to physical quadrature:

$$
\sum_K\int_K\nabla v\cdot\nabla u\,dx
\;\approx\;
\sum_{K,q}w_{Kq}\,
\nabla v|_K(x_{Kq})\cdot\nabla u|_K(x_{Kq}).
$$

The sum over $K$ has the two fixed local slots in this example. The package
handles mesh iteration and kernel dispatch. Weights contain physical measure. Reads retain cell identity, so two traces
at the same physical point can differ. Facet terms additionally need parent
cells, physical normals and facet weights.

The extension supplies its own domain, evaluation and derivative rules.
