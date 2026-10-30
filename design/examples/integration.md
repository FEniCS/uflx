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
one, two, three = Integer(1), Integer(2), Integer(3)
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
one = Integer(1)
e = domain_map(S, tensor((one - t, t)), T)
e_reversed = domain_map(S, tensor((t, one - t)), T)
F = tensor((Integer(0), X[0]))
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
e = domain_map(S, tensor((Integer(1) - t, t)), T)
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
zero, one = Integer(0), Integer(1)
F = tensor((zero, X[0]))
i, k = index("i"), index("k")

starts = tensor((tensor((zero, zero)), tensor((one, zero)), tensor((zero, one))))
ends = tensor((tensor((one, zero)), tensor((zero, one)), tensor((zero, zero))))
psi_edges = component_tensor(starts[i, k] + t * (ends[i, k] - starts[i, k]), (k,))
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
one, two, three = Integer(1), Integer(2), Integer(3)
psi = tensor((one + two * X[0] + X[1], X[0] + three * X[1]))
K = domain_map(T, psi, RD(2))
Y = coordinate(RD(2))
k = index("k")

area = integral(density(K), T)
centroid = component_tensor(integral(pull_back(Y[k], K) * density(K), T) / area, (k,))
```

The density is $5$, the area is $5/2$ and the centroid is $(2,4/3)$.
Since the map is affine, evaluating `psi` at the reference centroid
$(1/3,1/3)$ gives the same value through a different graph.

A decomposed domain permits the shorter vector-integral form:

```python
D = decomposed_domain((K,))
centroid_D = integral(coordinate(D), D) / integral(Integer(1), D)
```

`expand_integrals` splits the vector integrand into scalar components.

## Quadrature

An edge-midpoint rule integrates $f(x,y)=xy$ exactly on the unit triangle:

```python
T = UnitSimplex(2)
X = coordinate(T)
f = X[0] * X[1]
zero, half = Integer(0), Integer(1) / Integer(2)
sixth = Integer(1) / Integer(6)
i = index("i")

points = tensor((tensor((half, zero)), tensor((half, half)), tensor((zero, half))))
weights = tensor((sixth, sixth, sixth))
xs = domain_map(RD(0), points[i], T)
q = index_sum(weights[i] * integral(pull_back(f, xs), RD(0)), (i,))
```

The points form one family with free index `i`, so `f` occurs once in `q`.
The value is $1/24$, mathematically equal to `integral(f, T)`. Replacing the
integral by quadrature changes its graph and is a pass, with the rule supplied
as an argument.

### From quadrature to a loop

After lowering, evaluating the point integrals and reducing component tensors,
the expression is:

```python
index_sum(weights[i] * (points[i, 0] * points[i, 1]), (i,))
```

`expand_pull_backs` substitutes the points for `X`; the integral over `RD(0)`
returns its integrand. The compiler's `reduce_component_tensors` changes
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
| `RuntimeIndex(label, extent)` | An `AbstractIndexValue` of shape `()`, such as a local facet number; read with `gathered_index(value, ())`. |

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
cells = cell_domain(D, 0)
e, q = cell_index(T, 0), index("q")
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
