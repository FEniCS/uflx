# Surfaces and moving geometry

These examples use the core names from [language.md](../language.md) without
imports. Shared helper definitions are in the
[reference library](library.md).

## The Laplace-Beltrami operator

The same `gradient` writes Poisson's equation on a surface. Let `T` be a
topology of `ReferenceSimplex(2)` triangles, `x` a degree-1 coordinate
function of value shape `(3,)`, and `V` a real scalar Lagrange space on `T`.
The library supplies the topology, space and coordinate values:

```python
M = image(x, RD(3))

u, v, f = function(V, "u"), function(V, "v"), function(V, "f")
u_M, v_M, f_M = push_forward(u, M), push_forward(v, M), push_forward(f, M)

i = index("i")
grad_u = gradient(u_M, M)
grad_v = gradient(v_M, M)

a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), M), (v, u))
L = form(integral(f_M * v_M, M), (v,))
```

`M` has dimension 2 in `RD(3)`, so `grad_u` and `grad_v` are tangential
vectors with three components. Neither form needs an orientation.

On each reference cell, with $J$ the $3\times2$ coordinate Jacobian,
lowering gives the stiffness integrand

$$\widehat\nabla u^T (J^T J)^{-1}\widehat\nabla v\,
\sqrt{\det(J^T J)}.$$

The inverse metric comes from the pseudo-inverse in the gradient pull-back;
the square root is the surface density.

For second derivatives, set `P = gradient(coordinate(M), M)` and
`H = gradient(grad_u, M)`. `P` projects onto the tangent plane. `H` need not
be symmetric; the covariant Hessian is $P H P$.

```{dropdown} Legacy UFL
This is `inner(grad(u), grad(v)) * dx` on a mesh of topological dimension 2
in three dimensions. The domain dimensions determine the surface gradient.
```

## The flux through a surface

Flux through a closed surface needs an outward orientation; area does not.
Let `T` be a triangle topology with reference cell `K = ReferenceSimplex(2)`,
and `x` its coordinate function of value shape `(3,)`. The library provides
`T.orientation`, with signs $\pm1$ telling whether each cell's vertex order
agrees with the outward orientation.

```python
K = ReferenceSimplex(2)
M = image(x, RD(3))
cells = cell_domain(M, 0)
J = jacobian(cells)
r, l = index("r"), index("l")

def cross(a, b):
    return tensor((
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ))

nu = cross(component_tensor(J[r, 0], (r,)), component_tensor(J[r, 1], (r,)))

o = pull_back(T.orientation, cell_map(T, 0))

flux = index_sum(
    integral(o * index_sum(pull_back(coordinate(M), cells)[l] * nu[l], (l,)), K),
    (cell_index(T, 0),),
)

area = integral(Integer(1), M)
```

`nu` is the directed cofactor normal: its magnitude already includes the
surface density. Thus `flux` integrates on the reference cells without
another density factor. The signs `o` make every contribution outward, and
the position field's flux is three times the enclosed volume.

Consistently outward vertex ordering makes `o` unnecessary. Arbitrary ordering
without the signs instead lets contributions cancel incorrectly. `area` uses
$\sqrt{\det(J^T J)}$ and is independent of ordering. On a nonorientable
surface such as a Möbius strip, orientation is `None`: area remains defined,
but a global outward flux does not.

## A deformed domain

Move a mesh by changing its coordinate function. Functions on the topology
retain their degrees of freedom and push forward to the moved mesh.

For a two-triangle square, let the library supply the coordinate function `x`
with vertex values `((0,0), (1,0), (0,1), (1,1))`, a real vector displacement
space `W` of value shape `(2,)`, and a real scalar Lagrange space `V`, all on
`square`:

```python
K = ReferenceSimplex(2)
square = StaticTopology((K,), (((0, 1, 2), (1, 3, 2)),), ())
i = index("i")
s = function(W, "s")

D_s = image(x + s, RD(2))

volume = integral(Integer(1), D_s)

u, v = function(V, "u"), function(V, "v")
grad_u, grad_v = gradient(push_forward(u, D_s), D_s), gradient(push_forward(v, D_s), D_s)
a = form(integral(index_sum(grad_u[i] * grad_v[i], (i,)), D_s), (v, u))
```

For a cell `C = cell_domain(D_s, 0)`, the reference coordinate gradient is
$J+\widehat\nabla s$. Lowering therefore inserts:

| Quantity | Factor on the reference cell |
| -------- | ---------------------------- |
| Volume density | $\lvert\det(J+\widehat\nabla s)\rvert$ |
| Gradient pull-back | $(J+\widehat\nabla s)^{-T}$ |

Core derivatives hold domains fixed. The [shape-derivative extension](#a-shape-derivative)
first lowers the form, exposing `s` in these factors, then differentiates.
An argument as the direction gives a shape gradient.

A surface uses the same construction. For a triangle topology `T`, coordinate
function `x_T` of value shape `(3,)` and displacement space `W_T` on `T` of
the same value shape:

```python
w = function(W_T, "w")
M_w = image(x_T + w, RD(3))

area = integral(Integer(1), M_w)
```

The cell Jacobian is $3\times2$ and the density is $\sqrt{\det(J^T J)}$.

```{dropdown} Isogeometric analysis and shells
For extraction-based isogeometric analysis, store the NURBS map `F` as a
function on the topology of knot spans and use `image(F, RD(3))`. Multiple
patches use the combined topology; lowering inserts the geometric factors.

For a shell midsurface `image(phi_0, RD(3))`, the cell Jacobian is the
reference gradient of `phi_0`. Its Gram matrix is the first fundamental form,
and its density is the midsurface area element.
```

## A shape derivative

Shape differentiation must expose the geometry before differentiating.
Let `T` be a triangle topology, `X` a real vector Lagrange space of degree 1
and value shape `(2,)`, and `V` a real scalar Lagrange space, both on `T`.
The solver supplies coordinate values for `x`.
`shape_derivative` is the construction function of the
[extension](../language.md#shape-derivatives):

```python
x, h, v = function(X, "x"), function(X, "h"), function(X, "v")
D = image(x, RD(2))
grad_u = gradient(push_forward(function(V, "u"), D), D)
i = index("i")
energy = integral(index_sum(grad_u[i] * grad_u[i], (i,)), D)

d_energy = shape_derivative(energy, x, h)
g_energy = form(shape_derivative(energy, x, v), (v,))
```

`d_energy` is the derivative in the fixed direction `h`; `g_energy` binds an
argument and is the shape gradient. The extension runs in this order:

1. Lower the form so `x` appears in Jacobians and densities.
2. Replace `ShapeDerivative(g, x, h)` by `directional_derivative(g, x, x, h)`.
3. Apply derivatives and lower again.

Differentiating `energy` directly before step 1 gives `Zero`: `x` occurs only
inside the domain, which substitution does not enter. The extension needs
only the new node and this lowering rule; core differentiation does the rest.
