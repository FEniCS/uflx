# Geometry and integration

These examples use the core names from [language.md](../language.md) without
imports. Shared helper definitions are in the
[reference library](library.md).

## An affinely mapped triangle

Map the unit triangle to the triangle with vertices $(1,0)$, $(3,1)$ and
$(2,3)$:

$$\psi(X) = (1 + 2X_0 + X_1,\; X_0 + 3X_1).$$

Using the library's [UnitSimplex](library.md#the-unit-simplex):

```python
T = UnitSimplex(2)
X = coordinate(T)
one, two, three = integer(1), integer(2), integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))
K = domain_map(T, psi, RD(2))
```

The geometric quantities evaluate to:

| Expression | Value |
| ---------- | ----- |
| `jacobian(K)` | $\begin{pmatrix}2 & 1 \\ 1 & 3\end{pmatrix}$ |
| `determinant(jacobian(K))` | $5$ |
| `density(K)` | $5$ |

The physical area is therefore $5/2$. Swapping the last two vertices negates
the determinant but leaves the density and area unchanged.

## Circulation along an edge

Reversing an edge changes circulation's sign but leaves a scalar arc-length
integral unchanged. Take the field $F=(0,x)$ on the unit triangle and its
hypotenuse, traversed in both directions:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
one = integer(1)
e = domain_map(S, tensor((one - t, t)), T)
e_reversed = domain_map(S, tensor((t, one - t)), T)
F = tensor((integer(0), X[0]))
l = index("l")

def circulation(edge):
    return integral(index_sum(pull_back(F, edge)[l] * jacobian(edge)[l, 0], (l,)), S)
```

| Expression | Value after evaluation |
| ---------- | ---------------------- |
| `circulation(e)` | $1/2$ |
| `circulation(e_reversed)` | $-1/2$ |
| `integral(pull_back(X[0], e) * density(e), S)` | $\sqrt{2}/2$ |
| `integral(pull_back(X[0], e_reversed) * density(e_reversed), S)` | $\sqrt{2}/2$ |

The tangent is $(-1,1)$ on `e` and $(1,-1)$ on `e_reversed`. The circulation
integrands are respectively $1-t$ and $-t$, while both densities are
$\sqrt{2}$.

## Integrating over an edge

Flux uses a directed normal; arc length uses its magnitude. On the hypotenuse
of the unit triangle, the outward normal scaled by the edge length is
$(\psi'_1,-\psi'_0)$:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
e = domain_map(S, tensor((integer(1) - t, t)), T)
J = jacobian(e)

flux = integral(pull_back(X[0], e) * J[1, 0] - pull_back(X[1], e) * J[0, 0], S)
length = integral(density(e), S)
weighted_length = integral(pull_back(X[0] * X[1], e) * density(e), S)
```

The values are `flux` $=1$, `length` $=\sqrt{2}$ and `weighted_length`
$=\sqrt{2}/6$. The position field $X$ is tangent to the other two edges, so
all its outward flux passes through the hypotenuse. This agrees with the
divergence theorem: $\operatorname{div}X=2$ and the triangle has area $1/2$.

Here `density(e)` is $\sqrt{\det(J^T J)}$ with a $1\times1$ Gram matrix.
For a decomposed domain, `expand_integrals` inserts this density automatically.

## Integrating over the boundary of a triangle

Represent the three directed edges as one family, then sum their circulations:

```python
T, S = UnitSimplex(2), UnitSimplex(1)
X, t = coordinate(T), coordinate(S)[0]
zero, one = integer(0), integer(1)
F = tensor((zero, X[0]))
i, k = index("i"), index("k")

starts = tensor((tensor((zero, zero)), tensor((one, zero)), tensor((zero, one))))
ends = tensor((tensor((one, zero)), tensor((zero, one)), tensor((zero, zero))))
psi_edges = indexed_family(starts[i, k] + t * (ends[i, k] - starts[i, k]), (k,))
edges = domain_map(S, psi_edges, T)

integrand = index_sum(pull_back(F, edges)[k] * jacobian(edges)[k, 0], (k,))
result = index_sum(integral(integrand, S), (i,))
```

`result` evaluates to $1/2$, the circulation of $(0,x)$ around the triangle,
as Green's theorem predicts. Without the outer `index_sum`, the result has
free index `i`: one circulation per edge. All edges share the parameter domain
`S`; the order of `starts` and `ends` sets their directions.

## The centroid of a triangle

An integral is an expression, so a centroid can be written as a ratio of
integrals. For the triangle with vertices $(1,0)$, $(3,1)$ and $(2,3)$:

```python
T = UnitSimplex(2)
X = coordinate(T)
one, two, three = integer(1), integer(2), integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))
K = domain_map(T, psi, RD(2))
Y = coordinate(RD(2))
k = index("k")

area = integral(density(K), T)
centroid = indexed_family(integral(pull_back(Y[k], K) * density(K), T) / area, (k,))
```

The density is $5$, the area is $5/2$ and the centroid is $(2,4/3)$.
Since the map is affine, evaluating `psi` at the reference centroid
$(1/3,1/3)$ gives the same value through a different graph.

A decomposed domain permits the shorter vector-integral form:

```python
D = decomposed_domain((K,))
centroid_D = integral(coordinate(D), D) / integral(integer(1), D)
```

`expand_integrals` splits the vector integrand into scalar components.

## Integrating over point families

For explicit point families, use real, support-free `coordinates` of shape
`(n, 1)` and a support-free numerical array `a` of shape `(2,)`:

```python
Q = RD(1)
i = index("i")
p = indexed_family(domain_map(RD(0), coordinates[i], Q), (i,))
q = domain_map(RD(0), tensor((integer(3),)), Q)
points = decomposed_domain((p, q))
first_moment = integral(coordinate(Q)[0], points)
constant_integrals = integral(a[i], points)
```

With coordinate rows `(1,)` and `(2,)`, `first_moment` is 6 and
`constant_integrals` is `integer(3) * a[i]`, with free index `i` of extent 2.
The same index label is bound inside `p` and free in `a[i]`; neither domain
construction nor expansion captures it. The two pieces have ranks 1 and 0,
and every contribution retains that free `i`.

With two rows `(1,)` and `(1,)`, `first_moment` is 5: repeated points
contribute separately. For `n = 0`, `p` is empty, so only `q` contributes and
`first_moment` is 3. Using only the empty `p` gives zero, retaining the value
shape and free indices of any integrand. A family of shape `(n, m)` works
the same way with two fresh member indices, independent of an integrand's
free indices and of the indices bound inside that family.

## Located point samples

The mesh library supplies a two-cell interval mesh `D` with reference cell
`K = ReferenceSimplex(1)`. Its physical-cell parameterisations are
$X \mapsto X$ and $X \mapsto 1+X$. It supplies owner entries `(0, 1, 1)`, reference-position
rows `(1/2,)`, `(1/4,)`, `(1/4,)`, and weights `(2, 3, 5)` as support-free
families `owners`, `positions` and `weights` with no free indices.
`owners` has target extent 2; their shapes are `(3,)`, `(3, 1)` and `(3,)`.

```python
i = index("point")
f = coordinate(D)[0] ** integer(2)
reference_points = indexed_family(domain_map(RD(0), positions[i], K), (i,))
sampled = pull_back(pull_back(f, cell_domains(D, 0)[owners[i]]), reference_points[i])
value_i = integral(sampled, RD(0))
values = indexed_family(value_i, (i,))
weighted_total = index_sum(weights[i] * value_i, (i,))
```

The physical points are `(1/2, 5/4, 5/4)`, so `values` has entries
`(1/4, 25/16, 25/16)` and `weighted_total` is 13. The repeated point is
counted twice. Other free indices in `f` remain free; its scalar type and
argument spaces are retained. For zero points, the same constructions give
an empty sample array and zero total without evaluating a cell or point.

The [located-point contract](../language.md#a-family-of-points-in-a-mesh)
also applies after a compiler location pass, without a synthetic point
topology. Sampling `gradient(f, D)` gives the physical gradient at the
points; vector samples have one leading point axis.

## Quadrature

An edge-midpoint rule integrates $f(x,y)=xy$ exactly on the unit triangle:

```python
T = UnitSimplex(2)
X = coordinate(T)
f = X[0] * X[1]
zero, half = integer(0), integer(1) / integer(2)
sixth = integer(1) / integer(6)
i = index("i")

points = tensor((tensor((half, zero)), tensor((half, half)), tensor((zero, half))))
weights = tensor((sixth, sixth, sixth))
xs = indexed_family(domain_map(RD(0), points[i], T), (i,))
q = index_sum(weights[i] * integral(pull_back(f, xs[i]), RD(0)), (i,))
```

`xs` is a map family of shape `(3,)` with no free indices; `xs[i]` selects
one point with free index `i`. The expression `f` occurs once in `q`.
The value is $1/24$, mathematically equal to `integral(f, T)`. Replacing the
integral by quadrature changes its graph and is a pass, with the rule supplied
as an argument.

### From quadrature to a loop

After lowering, evaluating the point integrals and reducing indexed families,
the expression is:

```python
index_sum(weights[i] * (points[i, 0] * points[i, 1]), (i,))
```

Map-family lowering moves the selection outside the pull-back.
`expand_pull_backs` then substitutes the points for `X`; the integral over `RD(0)`
returns its integrand. The compiler's `reduce_indexed_families` changes
`points[i][0]` to `points[i, 0]`, and similarly for component 1. A generated
kernel could be:

```python
weights = (1 / 6, 1 / 6, 1 / 6)
points = ((1 / 2, 0.0), (1 / 2, 1 / 2), (0.0, 1 / 2))

def kernel() -> float:
    q = 0.0
    for i in range(3):
        q += weights[i] * (points[i][0] * points[i][1])
    return q
```

The `IndexSum` becomes the loop, and literals become tables. Runtime arrays
of shapes `(q, 2)` and `(q,)` use the same construction, with an
`AbstractExtent` `q` as the kernel's loop bound, as in the next example.

## Integrating over a cut cell

For a cell cut by $\phi<0$, the solver supplies quadrature points and weights
for the retained region. The number of points is known only at assembly.
These examples use the following library inputs, none of which belongs to the
core:

| Input | Meaning |
| ----- | ------- |
| `RuntimeExtent(label)` | An `AbstractExtent` known at run time. |
| `RuntimeArray(label, value_shape)` | An expression with no support holding quadrature data, coefficients or constants. |
| `RuntimeIndex(label, extent)` | A node extending `AbstractIndexValue`, exposed as an `IndexValue` of shape `()`; directly usable as a selector. |

With their lower-case construction functions, the cut-cell integral is:

```python
T = UnitSimplex(2)
X = coordinate(T)
f = X[0] * X[1]
i = index("i")
q = runtime_extent("q")
points = runtime_array("points", (q, 2))
weights = runtime_array("weights", (q,))

xs = domain_map(RD(0), points[i], T)
a = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The level set determines the solver's quadrature data; it does not appear in
the expression, since the core has no conditional. Point counts may differ
between cells. Representing all such rules in one topology-wide expression
requires the still-open [ragged indexing](../language.md#extents) convention.

## A quadrature-point coefficient

Material data can be stored at quadrature points without interpolating a field
across the cell. For a weighted mass form, supply one scalar coefficient value
per cell and quadrature point.

The library supplies a topology `T` of `ReferenceSimplex(2)` triangles, its
coordinate function `x` of value shape `(2,)`, a real scalar Lagrange space `V`
on `T`, and `runtime_array` for solver data. Use a three-point rule shared by
all cells; `points` and `weights` are reference-cell data, and `c` holds the
material values at the corresponding physical points:

```python
K = ReferenceSimplex(2)
D = image(x, RD(2))
e, q = index("e"), index("q")
cells = cell_domains(D, 0)[e]
points = runtime_array("points", (3, 2))
weights = runtime_array("weights", (3,))
c = runtime_array("c", (T.num_cells[0], 3))

u, v = function(V, "u"), function(V, "v")
sample = domain_map(RD(0), points[q], K)
u_q = pull_back(pull_back(push_forward(u, D), cells), sample)
v_q = pull_back(pull_back(push_forward(v, D), cells), sample)
rho_q = pull_back(density(cells), sample)

a = form(
    index_sum(integral(weights[q] * c[e, q] * u_q * v_q * rho_q, RD(0)), (e, q)),
    (v, u),
)
```

The inner pull-backs select a cell; the outer pull-backs sample on its reference points.
The density converts the reference weights to physical volume weights. The
sum over `q` forms the cell contribution, and the sum over `e` includes all
cells.

`c[e, q]` plays the role of a UFL quadrature-element coefficient. It has no
implied value between the points: changing their positions or ordering requires
matching data.
